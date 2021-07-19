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

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveTask;


public class Main {
	private static final String INPUT_FILENAME = "input";
	private static final int NUM_RUNS = 1;
	private static RawInput ri;

	public static void main(String[] args) {
		StopWatch sw = new StopWatch();
		sw.start();
//		System.out.close();

		for (int runNum = 0; runNum < NUM_RUNS; runNum++) {
			StopWatch lsw = new StopWatch();
			lsw.start();
			boolean isRunHessian = runNum % 2 == 0; // Hessian every other run

			AtomHandler.populateAtoms();
			InputHandler.processInput(INPUT_FILENAME);
			try {
				FileWriter fw1 = new FileWriter("errored-molecules.txt");
				fw1.write("");
				FileWriter fw2 = new FileWriter("dynamic-output.json");
				fw2.write("[");
				fw1.close();
				fw2.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			ri = InputHandler.ri;
			System.err.println(
					"MNDO Parameterization, updated 17 July 2021. " +
							ri.trainingSet + " training set (PM7)");
			System.err.println("input.json hash: " + ri.hash);

			// converts raw params array to NDDOParams classes and finds
			// params which need to be differentiated
			NDDOParams[] nddoParams = convertToNDDOParams(ri);

			// creates tasks to run in parallel and then runs them
			List<RecursiveTask> moleculeTasks = new ArrayList<>();
			for (RawMolecule request : ri.molecules) {
				if (request.isUsing) {
					moleculeTasks.add(new RecursiveTask() {
						@Override
						protected MoleculeRun compute() {
							MoleculeRun mr = new MoleculeRun(
									request,
									nddoParams,
									isRunHessian);
							mr.run();
							return mr;
						}
					});
				}
			}
			List<MoleculeRun> results = new ArrayList<>();
			for (ForkJoinTask task : ForkJoinTask.invokeAll(moleculeTasks))
				results.add((MoleculeRun) task.join()); // time intensive step

			// combined length of all differentiated params
			int paramLength = 0;
			for (int[] param : ri.neededParams) paramLength += param.length;

			// outputs output files for reference purposes, not read again in
			// this program
			List<MoleculeOutput> mos = new ArrayList<>(results.size());
			// processing results to change input.json for next run
			ParamOptimizer o = new ParamOptimizer();
			double[] ttGradient = new double[paramLength];
			double[][] ttHessian = isRunHessian ?
					new double[paramLength][paramLength] : null;
			for (MoleculeRun result : results) {
				// if it has not errored
				if (result.getRm().isUsing) {
					int[] moleculeATs = result.getRm().mats;
					int[][] moleculeNPs = result.getRm().mnps;
					boolean isDepad = true;

					initializePO(o, result, ri.neededParams, moleculeATs,
							moleculeNPs, isDepad);

					double[] g = ParamGradient.combine(
									result.getG().getTotalGradients(),
									ri.atomTypes, ri.neededParams,
									moleculeATs, moleculeNPs,
									isDepad);

					// ttGradient is sum of totalGradients across molecules
					for (int i = 0; i < g.length; i++) {
						ttGradient[i] += g[i];
					}

					if (isRunHessian) {
						double[][] h = result.getH().getHessian(
								ri.atomTypes, ri.neededParams);
						for (int i = 0; i < h.length; i++) {
							for (int j = 0; j < h[0].length; j++)
								ttHessian[i][j] += h[i][j];
						}
					}

					mos.add(OutputHandler.toMoleculeOutput(result));
				}
			}

			// optimizes params based on this run and gets new search direction
			DoubleMatrix newGradient = new DoubleMatrix(ttGradient);
			DoubleMatrix newHessian =
					isRunHessian ? new DoubleMatrix(ttHessian) :
							findMockHessian(newGradient,
									ri.params.lastHessian,
									ri.params.lastGradient,
									ri.params.lastDir,
									paramLength);
			double[] dir = o.optimize(newHessian, newGradient);

			// updating params and other input information
			ri.params.lastGradient = ttGradient;
			ri.params.lastHessian = ParamHessian.utify(
					newHessian.toArray2());
			ri.params.lastDir = dir;
			int n = 0;
			for (int atomI = 0; atomI < ri.neededParams.length; atomI++) {
				for (int paramI : ri.neededParams[atomI]) {
					ri.params.nddoParams[atomI][paramI] += dir[n];
					n++;
				}
			}

			// only place with actual file io
			try {
				Files.createDirectories(Path.of("outputs"));
				MoleculeOutput[] mosarray =
						mos.toArray(new MoleculeOutput[0]);
				OutputHandler.output(ri, mosarray, "outputs/run-"
						+ String.format("%04d", runNum) + "-output");
				InputHandler.updateInput(ri, INPUT_FILENAME);
			} catch (IOException e) {
				e.printStackTrace();
			}

			lsw.stop();
			System.err.println("\nRun " + runNum + " time taken: "
					+ lsw.getTime() + "\n\n---\n");
		}
		sw.stop();
		System.err.println("\nTotal time taken: " + sw.getTime());
	}

	private static void initializePO(ParamOptimizer o, MoleculeRun result,
									 int[][] neededParams,
									 int[] moleculeATs,
									 int[][] moleculeNP, boolean isDepad) {
		o.addData(new ReferenceData(result.getDatum()[0],
				result.getS().hf,
				ParamGradient.combine(
						result.getG().getHFDerivs(),
						ri.atomTypes, neededParams,
						moleculeATs, moleculeNP,
						isDepad),
				ReferenceData.HF_WEIGHT));
		if (result.getDatum()[1] != 0)
			o.addData(new ReferenceData(result.getDatum()[1],
					result.getS().dipole,
					ParamGradient.combine(
							result.getG().getDipoleDerivs(),
							ri.atomTypes, neededParams,
							moleculeATs, moleculeNP,
							isDepad),
					ReferenceData.DIPOLE_WEIGHT));
		if (result.getDatum()[2] != 0)
			o.addData(new ReferenceData(result.getDatum()[2],
					-result.getS().homo,
					ParamGradient.combine(
							result.getG().getIEDerivs(),
							ri.atomTypes, neededParams,
							moleculeATs, moleculeNP,
							isDepad),
					ReferenceData.IE_WEIGHT));
		if (result.isExpAvail()) o.addData(new ReferenceData(0,
				-result.getG().getE().getGeomGradient(),
				ParamGradient.combine(
						result.getG().getGeomDerivs(),
						ri.atomTypes, neededParams,
						moleculeATs, moleculeNP,
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

	private static NDDOParams[] convertToNDDOParams(RawInput ri) {
		NDDOParams[] nddoParams = null;
		switch (ri.model) {
			case "mndo":
				nddoParams = new MNDOParams[ri.params.nddoParams.length];
				for (int i = 0; i < ri.params.nddoParams.length; i++)
					nddoParams[i] = new MNDOParams(ri.params.nddoParams[i]);
				break;
			case "am1":
				nddoParams = new AM1Params[ri.params.nddoParams.length];
				for (int i = 0; i < ri.params.nddoParams.length; i++)
					nddoParams[i] = new AM1Params(ri.params.nddoParams[i]);
				break;
		}

		assert nddoParams != null;
		return nddoParams;
	}
}