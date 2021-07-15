import nddoparam.NDDOParams;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamHessian;
import optimize.ParamOptimizer;
import optimize.ReferenceData;
import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import runcycle.MoleculeRun;
import runcycle.MoleculeRunR;
import runcycle.MoleculeRunU;
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

	public static void main(String[] args) {
		StopWatch sw = new StopWatch();
		sw.start();
//        System.out.close();

		for (int numRuns = 0; numRuns < 1; numRuns++) {
			boolean runHessian = numRuns % 2 == 0; // Hessian every other run

			AtomHandler.populateAtoms();
			InputHandler.processInput(INPUT_FILENAME);

			RawInput ri = InputHandler.ri;
			System.out.println(
					"MNDO Parameterization, updated 13 July. " +
							ri.trainingSet +
							" training set (PM7)");

			NDDOParams[] nddoParams =
					new MNDOParams[ri.params.nddoParams.length];
			// TODO change the following line if AM1
			for (int i = 0; i < ri.params.nddoParams.length; i++)
				nddoParams[i] = new MNDOParams(ri.params.nddoParams[i]);

			int[][] neededParams = new int[ri.atomTypes.length][];
			int w = 0;
			for (int atomType : ri.atomTypes) {
				if (atomType == 1)
					neededParams[w] = MNDOParams.T1ParamNums;
				else neededParams[w] = MNDOParams.T2ParamNums;
				w++;
			}

			try {
				List<RawMolecule> requests =
						new ArrayList<>(Arrays.asList(ri.molecules));

				int cores = Runtime.getRuntime().availableProcessors();
				int remainingNonParallel = 5;
				int maxParallel = remainingNonParallel < requests.size() ?
						requests.size() - remainingNonParallel : 1;
				List<RawMolecule> parallelRequests =
						requests.subList(0, maxParallel);
				ForkJoinPool threadPool = new ForkJoinPool(cores);

				List<MoleculeRun> results = threadPool
						.submit(() -> parallelRequests.parallelStream()
								.map(request -> {
									MoleculeRun result = request.restricted ?
											new MoleculeRunR(request,
													nddoParams,
													ri.atomTypes, runHessian) :
											new MoleculeRunU(request,
													nddoParams,
													ri.atomTypes, runHessian);
									return result;
								})).get().collect(Collectors.toList());

				for (RawMolecule request : requests
						.subList(maxParallel, requests.size())) {
					MoleculeRun result = request.restricted ?
							new MoleculeRunR(request, nddoParams, ri.atomTypes,
									runHessian) :
							new MoleculeRunU(request, nddoParams,
									ri.atomTypes,
									runHessian);
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
				if (runHessian) ttHessian =
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
					if (runHessian) {
						double[][] h = result.getH().getHessianUnpadded();
						for (int i = 0; i < h.length; i++) {
							for (int j = 0; j < h[0].length; j++)
								ttHessian[i][j] += h[i][j];
						}
					}
				}
				DoubleMatrix newGradient = new DoubleMatrix(ttGradient);
				DoubleMatrix B = runHessian ? new DoubleMatrix(ttHessian) :
						findMockHessian(newGradient,
								ri.params.lastHessian,
								ri.params.lastGradient, ri.params.lastDir,
								paramLength);
				double[] dir = o.optimize(B, newGradient);

				ri.params.lastGradient  = ttGradient;
				ri.params.lastHessian = ParamHessian.getHessianUT(ttHessian);
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

			} catch (Exception e) {
				e.printStackTrace();
			}
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