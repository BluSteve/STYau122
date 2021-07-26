import nddoparam.NDDOParams;
import nddoparam.am1.AM1Params;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamGradient;
import nddoparam.param.ParamHessian;
import optimize.ParamOptimizer;
import optimize.ReferenceData;
import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import runcycle.MoleculeRan;
import runcycle.MoleculeResult;
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
import java.util.*;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.RecursiveTask;


public class Main {
	private static final String INPUT_FILENAME = "input";
	private static final int NUM_RUNS = 1;
	private static final boolean isImportLastRun = true;
	private static final int MAX_RETRIES = 2;
	private static RawInput ri;
	private static MoleculeOutput[] ranMolecules;
	private static boolean isRetryFailed = false;

	public static void main(String[] args) {
		StopWatch sw = new StopWatch();
		sw.start();
		Scanner s = new Scanner(System.in);
		System.err.print("Run in verbose mode? (y/N) ");
		if (!s.next().equals("y"))
			System.out.close();

		if (isImportLastRun) {
			ranMolecules =
					OutputHandler.importMoleculeOutputs("dynamic-output");
			if (ranMolecules != null) {
				System.err.print(ranMolecules.length +
						" molecules from last run found, would you like to " +
						"import them? (Y/n) ");
				if (s.next().equals("n")) {
					System.err.println("Last ran molecules import cancelled!");
					ranMolecules = null;
				}
				else {
					System.err.println(ranMolecules.length +
							" molecules successfully imported from last" +
							" incomplete run!");
				}
			}
			else System.err.println("Last ran molecules import failed!");
			System.err.println();
		}

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
			System.err.println("Run number: " + runNum);

			// converts raw params array to NDDOParams classes and finds
			// params which need to be differentiated
			NDDOParams[] nddoParams = convertToNDDOParams(ri);
			// combined length of all differentiated params
			int paramLength = 0;
			for (int[] param : ri.neededParams) paramLength += param.length;

			// creates tasks to run in parallel and then runs them
			int failedCount = 0;
			List<RecursiveTask> moleculeTasks = new ArrayList<>();
			for (RawMolecule request : ri.molecules) {
				if (request.isUsing) {
					boolean isDoneAlready = false;
					if (ranMolecules != null) {
						for (MoleculeOutput mo : ranMolecules) {
							if (mo.rawMolecule.index == request.index) {
								isDoneAlready = true;
								break;
							}
						}
					}
					if (!isDoneAlready) {
						moleculeTasks.add(new RecursiveTask() {
							@Override
							protected MoleculeRun compute() {
								MoleculeRun mr = new MoleculeRun(
										request, nddoParams, isRunHessian);
								mr.run();
								return mr;
							}
						});
					}
				}
				else failedCount++;
			}

			System.err.println(
					failedCount + " previously failed molecules skipped");

			Collections.shuffle(moleculeTasks, new Random(123));
//			moleculeTasks = moleculeTasks.subList(0,10);
			List<MoleculeResult> results = new ArrayList<>();
			for (ForkJoinTask task : ForkJoinTask.invokeAll(moleculeTasks)) {
				results.add((MoleculeRun) task.join());
			}
			if (ranMolecules != null) {
				for (MoleculeOutput mr : ranMolecules) {
					results.add(new MoleculeRan(mr));
				}
			}

			// outputs output files for reference purposes, not read again in
			// this program
			List<MoleculeOutput> mos = new ArrayList<>(results.size());

			// processing results to change input.json for next run
			ParamOptimizer o = new ParamOptimizer();
			double[] ttGradient = new double[paramLength];
			double[][] ttHessian = isRunHessian ?
					new double[paramLength][paramLength] : null;

			if (isRetryFailed) {
				// retries failed molecules until they succeed
				failedCount = 0;
				for (MoleculeResult result : results) {
					if (!result.getRm().isUsing) failedCount++;
				}
				int retries = 0;
				while (failedCount > 0 && retries < MAX_RETRIES) {
					System.err.println(
							"Rerunning " + failedCount +
									" failed molecules for the " + retries +
									"th time.");
					List<RecursiveAction> moleculeReruns = new ArrayList<>();
					for (MoleculeResult result : results) {
						// if isUsing is false here then that means the
						// molecule failed during this run, not a previous one
						if (result instanceof MoleculeRun &&
								!result.getRm().isUsing) {
							result.getRm().isUsing = true;
							moleculeReruns.add(new RecursiveAction() {
								@Override
								protected void compute() {
									((MoleculeRun) result).run();
								}
							});

						}
					}
					ForkJoinTask.invokeAll(moleculeReruns);
					failedCount = 0;
					for (MoleculeResult result : results) {
						if (!result.getRm().isUsing) failedCount++;
					}
					retries++;
				}
			}

			for (MoleculeResult result : results) {
				// if it has not errored
				if (result.getRm().isUsing) {
					int[] moleculeATs = result.getRm().mats;
					int[][] moleculeNPs = result.getRm().mnps;
					boolean isDepad = true;

					addToPO(o, result, ri.neededParams, moleculeATs,
							moleculeNPs, isDepad);

					double[] g = ParamGradient.combine(
							result.getTotalGradients(),
							ri.atomTypes, ri.neededParams,
							moleculeATs, moleculeNPs,
							isDepad);

					// ttGradient is sum of totalGradients across molecules
					for (int i = 0; i < g.length; i++) {
						ttGradient[i] += g[i];
					}

					if (isRunHessian) {
						double[][] h =
								ParamHessian.padHessian(result.getHessian(),
										result.getRm().mats,
										ri.atomTypes, ri.neededParams);
						for (int i = 0; i < h.length; i++) {
							for (int j = 0; j < h[0].length; j++)
								ttHessian[i][j] += h[i][j];
						}
					}

					mos.add(OutputHandler.toMoleculeOutput(result, isRunHessian));
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

	private static void addToPO(ParamOptimizer o, MoleculeResult result,
								int[][] neededParams,
								int[] moleculeATs,
								int[][] moleculeNP, boolean isDepad) {
		o.addData(new ReferenceData(result.getDatum()[0],
				result.getHF(),
				ParamGradient.combine(
						result.getHFDerivs(),
						ri.atomTypes, neededParams,
						moleculeATs, moleculeNP,
						isDepad),
				ReferenceData.HF_WEIGHT));
		if (result.getDatum()[1] != 0)
			o.addData(new ReferenceData(result.getDatum()[1],
					result.getDipole(),
					ParamGradient.combine(
							result.getDipoleDerivs(),
							ri.atomTypes, neededParams,
							moleculeATs, moleculeNP,
							isDepad),
					ReferenceData.DIPOLE_WEIGHT));
		if (result.getDatum()[2] != 0)
			o.addData(new ReferenceData(result.getDatum()[2],
					result.getIE(),
					ParamGradient.combine(
							result.getIEDerivs(),
							ri.atomTypes, neededParams,
							moleculeATs, moleculeNP,
							isDepad),
					ReferenceData.IE_WEIGHT));
		if (result.isExpAvail()) o.addData(new ReferenceData(0,
				result.getGeomGradient(),
				ParamGradient.combine(
						result.getGeomDerivs(),
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