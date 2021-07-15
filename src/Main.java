import nddoparam.NDDOParams;
import nddoparam.am1.AM1Params;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamHessian;
import optimize.ParamOptimizer;
import optimize.ReferenceData;
import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import runcycle.MoleculeRun;
import runcycle.input.InputHandler;
import runcycle.input.RawInput;
import runcycle.input.RawMolecule;
import runcycle.output.MoleculeOutput;
import runcycle.output.OutputHandler;
import scf.AtomHandler;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;


public class Main {

	private static final String INPUT_FILENAME = "input.json";
	private static final String OUTPUT_FILENAME = "output.json";
	private static RawInput ri;

	public static void main(String[] args) {
		StopWatch sw = new StopWatch();
		sw.start();
//        System.out.close();

		for (int numRuns = 0; numRuns < 1; numRuns++) {
			boolean isRunHessian = numRuns % 2 == 0; // Hessian every other run

			AtomHandler.populateAtoms();
			InputHandler.processInput(INPUT_FILENAME);

			ri = InputHandler.ri;
			System.out.println(
					"MNDO Parameterization, updated 13 July. " +
							ri.trainingSet +
							" training set (PM7)");

			NDDOParams[] nddoParams = null;

			switch (ri.model) {
				case "mndo":
					nddoParams =
							new MNDOParams[ri.params.nddoParams.length];
					for (int i = 0; i < ri.params.nddoParams.length; i++)
						nddoParams[i] =
								new MNDOParams(ri.params.nddoParams[i]);
					break;
				case "am1":
					nddoParams =
							new AM1Params[ri.params.nddoParams.length];
					for (int i = 0; i < ri.params.nddoParams.length; i++)
						nddoParams[i] = new AM1Params(ri.params.nddoParams[i]);
					break;
			}


			int[][] neededParams = new int[ri.atomTypes.length][];
			int w = 0;
			for (int atomType : ri.atomTypes) {
				switch (ri.model) {
					case "mndo":
						if (atomType == 1)
							neededParams[w] = MNDOParams.T1ParamNums;
						else neededParams[w] = MNDOParams.T2ParamNums;
						break;
					case "am1":
						if (atomType == 1)
							neededParams[w] = AM1Params.T1ParamNums;
						else neededParams[w] = AM1Params.T2ParamNums;
						break;
				}

				w++;
			}

			List<RawMolecule> requests =
					new ArrayList<>(Arrays.asList(ri.molecules));

			int cores = Runtime.getRuntime().availableProcessors();
			int remainingNonParallel = 0;
			int maxParallel = remainingNonParallel < requests.size() ?
					requests.size() - remainingNonParallel : 1;
			List<RawMolecule> parallelRequests =
					requests.subList(0, maxParallel);
			ForkJoinPool threadPool = new ForkJoinPool(cores);

			List<MoleculeRun> results = null;
			try {
				NDDOParams[] finalNddoParams = nddoParams;
				results = threadPool
						.submit(() -> parallelRequests.parallelStream()
								.map(request -> new MoleculeRun(
										request,
										finalNddoParams,
										ri.atomTypes,
										isRunHessian)))
						.get()
						.collect(Collectors.toList());
			}
			catch (Exception e) {
				e.printStackTrace();
			}

			for (RawMolecule rawMolecule : requests
					.subList(maxParallel, requests.size())) {
				MoleculeRun result =
						new MoleculeRun(rawMolecule, nddoParams,
								ri.atomTypes, isRunHessian);
				results.add(result);
			}

			MoleculeOutput[] mos = new MoleculeOutput[results.size()];
			for (int i = 0; i < results.size(); i++) {
				mos[i] = OutputHandler.toMoleculeOutput(results.get(i));
			}

			ParamOptimizer o = new ParamOptimizer();
			int paramLength =
					results.get(0).getG().combine(results.get(0).getG()
							.depad(results.get(0).getG()
									.getTotalGradients())).length;
			double[] ttGradient = new double[paramLength];
			double[][] ttHessian = null;
			if (isRunHessian) ttHessian =
					new double[results.get(0).getH()
							.getHessianUnpadded().length]
							[results.get(0).getH()
							.getHessianUnpadded()[0].length];
			for (MoleculeRun result : results) {
				o.addData(new ReferenceData(result.getDatum()[0],
						result.getS().hf,
						result.getG().combine(
								result.getG().depad(result.getG()
										.getHFDerivs())),
						ReferenceData.HF_WEIGHT));
				if (result.getDatum()[1] != 0)
					o.addData(new ReferenceData(result.getDatum()[1],
							result.getS().dipole,
							result.getG().combine(
									result.getG().depad(result.getG()
											.getDipoleDerivs())),
							ReferenceData.DIPOLE_WEIGHT));
				if (result.getDatum()[2] != 0)
					o.addData(new ReferenceData(result.getDatum()[2],
							-result.getS().homo,
							result.getG()
									.combine(result.getG()
											.depad(result.getG()
													.getIEDerivs())),
							ReferenceData.IE_WEIGHT));
				if (result.isExpAvail()) o.addData(new ReferenceData(0,
						-result.getG().getE().geomGradient,
						result.getG()
								.combine(result.getG().depad(result.getG()
										.getGeomDerivs())),
						ReferenceData.GEOM_WEIGHT));

				double[] g =
						result.getG().combine(result.getG()
								.depad(result.getG().getTotalGradients()));
				for (int i = 0; i < g.length; i++) {
					ttGradient[i] += g[i];
				}
				if (isRunHessian) {
					double[][] h = result.getH().getHessianUnpadded();
					for (int i = 0; i < h.length; i++) {
						for (int j = 0; j < h[0].length; j++)
							ttHessian[i][j] += h[i][j];
					}
				}
			}
			DoubleMatrix newGradient = new DoubleMatrix(ttGradient);
			DoubleMatrix newHessian =
					isRunHessian ? new DoubleMatrix(ttHessian) :
							findMockHessian(newGradient,
									ri.params.lastHessian,
									ri.params.lastGradient,
									ri.params.lastDir,
									paramLength);
			double[] dir = o.optimize(newHessian, newGradient);

			ri.params.lastGradient = ttGradient;
			ri.params.lastHessian = ParamHessian.getHessianUT(
					newHessian.toArray2());
			ri.params.lastDir = dir;

			int n = 0;
			for (int atomI = 0; atomI < neededParams.length; atomI++) {
				for (int paramI : neededParams[atomI]) {
					ri.params.nddoParams[atomI][paramI] += dir[n];
					n++;
				}
			}

			OutputHandler.output(mos, OUTPUT_FILENAME);
			InputHandler.updateInput(ri, INPUT_FILENAME);

		}

		sw.stop();
		System.out.println("\nTime taken: " + sw.getTime());
	}

	private static DoubleMatrix findMockHessian(DoubleMatrix newGradient,
												double[] oldHessian,
												double[] oldGradient,
												double[] oldDir, int size) {
		DoubleMatrix s = new DoubleMatrix(oldDir);
		DoubleMatrix hessian = new DoubleMatrix(size, size);
		int count = 1;
		int index = 0;
		while (index < ((size + 1) * size) / 2) {
			for (int i = count - 1; i < size; i++) {
				hessian.put(count - 1, i, oldHessian[index]);
				hessian.put(i, count - 1, oldHessian[index]);
				index++;
			}
			count++;
		}
		DoubleMatrix y = newGradient.sub(new DoubleMatrix(oldGradient));

		double b = y.transpose().mmul(s).get(0);
		DoubleMatrix A = y.mmul(y.transpose()).mmul(1 / b);
		double a = s.transpose().mmul(hessian).mmul(s).get(0);
		DoubleMatrix C = hessian.mmul(s).mmul(s.transpose()).mmul
				(hessian.transpose()).mmul(1 / a);

		return hessian.add(A).sub(C);
	}

}