import nddoparam.NDDOParams;
import nddoparam.am1.AM1Params;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamGradient;
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
import java.util.stream.Collectors;


public class Main {

	private static final String INPUT_FILENAME = "input.json";
	private static final int NUM_RUNS = 1;
	private static RawInput ri;

	public static void main(String[] args) {
		StopWatch sw = new StopWatch();
		sw.start();
		System.out.close();

		for (int runNum = 0; runNum < NUM_RUNS; runNum++) {
			StopWatch lsw = new StopWatch();
			lsw.start();
			boolean isRunHessian = runNum % 2 == 0; // Hessian every other run

			AtomHandler.populateAtoms();
			InputHandler.processInput(INPUT_FILENAME);

			ri = InputHandler.ri;
			System.err.println(
					"MNDO Parameterization, updated 16 July. " +
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

			List<MoleculeRun> results = null;
			try {
				NDDOParams[] finalNddoParams = nddoParams;
				results = requests.parallelStream()
						.map(request -> new MoleculeRun(
								request,
								finalNddoParams,
								ri.atomTypes,
								isRunHessian)).collect(Collectors.toList());
			} catch (Exception e) {
				e.printStackTrace();
			}

			ParamOptimizer o = new ParamOptimizer();
			// combined length of all differentiated params
			int paramLength = 0;
			for (int[] param : neededParams) paramLength += param.length;

			double[] ttGradient = new double[paramLength];
			double[][] ttHessian = null;
			// todo pad hessian
			if (isRunHessian) ttHessian = new double[paramLength][paramLength];

			for (MoleculeRun result : results) {
				int[] moleculeUZ = result.getS().getUniqueZs();
				int[][] moleculeNP = result.getS().getNeededParams();
				boolean isDepad = true;

				initializePO(neededParams, o, result, moleculeUZ, moleculeNP,
						isDepad);

				double[] g =
						ParamGradient.combine(
								result.getG().getTotalGradients(),
								ri.atomTypes, neededParams,
								moleculeUZ, moleculeNP,
								isDepad);

				// ttGradient is sum of totalGradients across molecules
				for (int i = 0; i < g.length; i++) {
					ttGradient[i] += g[i];
				}
				if (isRunHessian) {
					double[][] h =
							result.getH().getHessianUnpadded(ri.atomTypes,
									neededParams, paramLength);
					for (int i = 0; i < h.length; i++) {
						for (int j = 0; j < h[0].length; j++)
							ttHessian[i][j] += h[i][j];
					}
					System.err.println(Arrays.deepToString(h));

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


			MoleculeOutput[] mos = new MoleculeOutput[results.size()];
			for (int i = 0; i < results.size(); i++) {
				mos[i] = OutputHandler.toMoleculeOutput(results.get(i));
			}
			OutputHandler.output(mos, "output.json");
			OutputHandler.output(mos,
					"outputs/run-" + String.format("%04d", runNum) +
							"-output.json");

			InputHandler.updateInput(ri, INPUT_FILENAME);

			System.err.println(
					"\nRun " + runNum + " time taken: " + lsw.getTime() +
							"\n\n---\n");
		}

		sw.stop();
		System.err.println("\nTotal time taken: " + sw.getTime());
	}

	private static void initializePO(int[][] neededParams, ParamOptimizer o,
									 MoleculeRun result, int[] moleculeUZ,
									 int[][] moleculeNP, boolean isDepad) {
		o.addData(new ReferenceData(result.getDatum()[0],
				result.getS().hf,
				ParamGradient.combine(
						result.getG().getHFDerivs(),
						ri.atomTypes, neededParams,
						moleculeUZ, moleculeNP,
						isDepad),
				ReferenceData.HF_WEIGHT));
		if (result.getDatum()[1] != 0)
			o.addData(new ReferenceData(result.getDatum()[1],
					result.getS().dipole,
					ParamGradient.combine(
							result.getG().getDipoleDerivs(),
							ri.atomTypes, neededParams,
							moleculeUZ, moleculeNP,
							isDepad),
					ReferenceData.DIPOLE_WEIGHT));
		if (result.getDatum()[2] != 0)
			o.addData(new ReferenceData(result.getDatum()[2],
					-result.getS().homo,
					ParamGradient.combine(
							result.getG().getIEDerivs(),
							ri.atomTypes, neededParams,
							moleculeUZ, moleculeNP,
							isDepad),
					ReferenceData.IE_WEIGHT));
		if (result.isExpAvail()) o.addData(new ReferenceData(0,
				-result.getG().getE().geomGradient,
				ParamGradient.combine(
						result.getG().getGeomDerivs(),
						ri.atomTypes, neededParams,
						moleculeUZ, moleculeNP,
						isDepad),
				ReferenceData.GEOM_WEIGHT));
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