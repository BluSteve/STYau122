import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import optimize.ParamOptimizer;
import optimize.ReferenceData;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.MoleculeRan;
import runcycle.MoleculeResult;
import runcycle.MoleculeRun;
import runcycle.input.InputHandler;
import runcycle.input.RawInput;
import runcycle.input.RawMolecule;
import runcycle.output.MoleculeOutput;
import runcycle.output.OutputHandler;
import scf.AtomHandler;
import tools.Utils;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;


public class Main {
	private static final String INPUT_FILENAME = "input";
	private static final int NUM_RUNS = 270;
	private static final int MAX_MOLECULES = 1000;
	private static final boolean isImportLastRun = true;
	private static final Logger logger = LogManager.getLogger();
	private static final ScheduledExecutorService progressBar =
			Executors.newScheduledThreadPool(1);
	private static RawInput ri;
	private static MoleculeOutput[] ranMolecules;

	public static void main(String[] args) throws IOException {
		System.out.println(logger.getLevel());

		if (isImportLastRun) {
			ranMolecules = OutputHandler.importMoleculeOutputs("dynamic-output");
			if (ranMolecules != null) {
				Scanner s = new Scanner(System.in);
				System.out.print(ranMolecules.length + " molecules from last run found, would you like to " +
						"import them? (Y/n) ");
				if (s.next().equals("n")) {
					logger.info("Last ran molecules import cancelled!");
					ranMolecules = null;
				}
				else {
					logger.info("{} molecules successfully imported from last incomplete run!", ranMolecules.length);
				}
				s.close();
			}
			else logger.warn("Last ran molecules import failed!");
		}

		logger.info("NUM_RUNS = " + NUM_RUNS);

		boolean[] isDones = new boolean[MAX_MOLECULES];
		if (logger.isInfoEnabled()) {
			Runnable mLeft = () -> {
				List<String> done = new ArrayList<>();
				List<String> left = new ArrayList<>();

				for (RawMolecule rm : ri.molecules) {
					if (isDones[rm.index]) done.add(rm.debugName());
					else left.add(rm.debugName());
				}

				done.sort(String::compareTo);
				left.sort(String::compareTo);

				logger.info("Molecules done ({}): {}", done.size(), done);
				logger.info("Molecules left ({}): {}", left.size(), left);
			};

			int wait = 10;
			progressBar.scheduleAtFixedRate(mLeft, wait, wait,
					TimeUnit.SECONDS);
		}

		StopWatch sw = new StopWatch();
		sw.start();

		for (int runNum = 0; runNum < NUM_RUNS; runNum++) {
			StopWatch lsw = new StopWatch();
			lsw.start();


			boolean isRunHessian = runNum % 2 == 0; // Hessian every other run

			AtomHandler.populateAtoms();
			ri = InputHandler.processInput(INPUT_FILENAME);

			Arrays.fill(isDones, false);

			FileWriter fw2 = new FileWriter("dynamic-output.json");
			fw2.write("[");
			fw2.close();

			if (runNum == 0) {
				logger.info("MNDO Parameterization, updated 21 Nov 2021. {} training set (PM7)", ri.trainingSet);
			}
			logger.info("Run number: {}, hessian: {}, {}.json hash: {}\n", runNum, isRunHessian, INPUT_FILENAME,
					ri.hash);


			// this array is used as a Map<Integer, NDDOParams>
			NDDOParams[] npMap = Utils.getNpMap(ri);

			// combined length of all differentiated params
			int paramLength = 0;
			for (int[] param : ri.neededParams) paramLength += param.length;

			// create tasks to run in parallel and then runs them excludes ran molecules from previous dynamic output
			List<RecursiveTask<MoleculeRun>> moleculeTasks = new ArrayList<>();
			for (RawMolecule rm : ri.molecules) {
				boolean isDoneAlready = false;
				if (ranMolecules != null) {
					for (MoleculeOutput mo : ranMolecules) {
						if (mo.rawMolecule.index == rm.index) {
							isDoneAlready = true;
							isDones[rm.index] = true;
							break;
						}
					}
				}
				if (!isDoneAlready) {
					moleculeTasks.add(new RecursiveTask<>() {
						@Override
						protected MoleculeRun compute() {
							NDDOAtom[] atoms = Utils.toNDDOAtoms(ri.model, rm.atoms, npMap);
							NDDOAtom[] expGeom = rm.expGeom == null ? null :
									Utils.toNDDOAtoms(ri.model, rm.expGeom, npMap);

							MoleculeRun mr = new MoleculeRun(rm, atoms, expGeom, isRunHessian);
							mr.run();

							isDones[rm.index] = true;
							return mr;
						}
					});
				}
			}

			List<MoleculeResult> results = new ArrayList<>();

			// adds previously ran molecules into the final results list
			if (ranMolecules != null) {
				for (MoleculeOutput mr : ranMolecules) {
					results.add(new MoleculeRan(mr));
				}
			}

			// shuffles run requests and runs them
			Collections.shuffle(moleculeTasks, new Random(123));
			for (ForkJoinTask<MoleculeRun> task : ForkJoinTask.invokeAll(moleculeTasks)) {
				results.add(task.join());
			}

			// outputs output files for reference purposes, not read again in this program
			List<MoleculeOutput> mos = new ArrayList<>(results.size());

			// processing results to change input.json for next run
			ParamOptimizer o = new ParamOptimizer();
			double[] ttGradient = new double[paramLength];
			double[][] ttHessian = isRunHessian ? new double[paramLength][paramLength] : null;

			for (MoleculeResult result : results) {
				int[] moleculeATs = result.getRm().mats;
				int[][] moleculeNPs = result.getRm().mnps;
				boolean isDepad = true;

				addToPO(o, result, ri.neededParams, moleculeATs, moleculeNPs, isDepad);

				double[] g = ParamGradient.combine(result.getTotalGradients(), ri.atomTypes, ri.neededParams,
						moleculeATs, moleculeNPs, isDepad);

				// ttGradient is sum of totalGradients across molecules
				for (int i = 0; i < g.length; i++) {
					ttGradient[i] += g[i];
				}

				if (isRunHessian) {
					double[][] h = ParamHessian.padHessian(result.getHessian(), result.getRm().mats, ri.atomTypes,
							ri.neededParams);

					for (int i = 0; i < h.length; i++) {
						for (int j = 0; j < h[0].length; j++)
							ttHessian[i][j] += h[i][j];
					}
				}

				mos.add(OutputHandler.toMoleculeOutput(result, isRunHessian));
			}

			mos.sort(Comparator.comparingInt(x -> x.rawMolecule.index));

			// optimizes params based on this run and gets new search direction
			SimpleMatrix newGradient = new SimpleMatrix(ttGradient);
			SimpleMatrix newHessian = isRunHessian ?
					new SimpleMatrix(ttHessian) :
					findMockHessian(newGradient, ri.params.lastHessian, ri.params.lastGradient, ri.params.lastDir,
							paramLength);

			double[] dir = o.optimize(newHessian, newGradient);

			// updating params and other input information
			int n = 0;
			for (int atomI = 0; atomI < ri.neededParams.length; atomI++) {
				for (int paramI : ri.neededParams[atomI]) {
					ri.params.nddoParams[atomI][paramI] += dir[n];
					n++;
				}
			}
			ri.params.lastGradient = ttGradient;
			ri.params.lastHessian = ParamHessian.utify(newHessian.toArray2());
			ri.params.lastDir = dir;

			lsw.stop();
			logger.info("Run {} time taken: {}", runNum, lsw.getTime());

			Files.createDirectories(Path.of("outputs"));
			MoleculeOutput[] mosarray = mos.toArray(new MoleculeOutput[0]);
			OutputHandler.output(ri, mosarray, "outputs/run-" + String.format("%04d", runNum) + "-output",
					lsw.getTime());
			InputHandler.outputInput(ri, INPUT_FILENAME);
		}
		sw.stop();
		logger.info("Total time taken: {}", sw.getTime());

		// stop progress bar and terminate program
		progressBar.shutdownNow();
		System.exit(0);
	}

	private static void addToPO(ParamOptimizer o, MoleculeResult result, int[][] neededParams, int[] moleculeATs,
								int[][] moleculeNP, boolean isDepad) {
		o.addData(new ReferenceData(result.getDatum()[0], result.getHF(),
				ParamGradient.combine(result.getHFDerivs(), ri.atomTypes, neededParams, moleculeATs, moleculeNP,
						isDepad), ReferenceData.HF_WEIGHT));

		if (result.getDatum()[1] != 0) o.addData(new ReferenceData(result.getDatum()[1], result.getDipole(),
				ParamGradient.combine(result.getDipoleDerivs(), ri.atomTypes, neededParams, moleculeATs, moleculeNP,
						isDepad), ReferenceData.DIPOLE_WEIGHT));

		if (result.getDatum()[2] != 0) o.addData(new ReferenceData(result.getDatum()[2], result.getIE(),
				ParamGradient.combine(result.getIEDerivs(), ri.atomTypes, neededParams, moleculeATs, moleculeNP,
						isDepad), ReferenceData.IE_WEIGHT));

		if (result.isExpAvail()) o.addData(new ReferenceData(0, result.getGeomGradient(),
				ParamGradient.combine(result.getGeomDerivs(), ri.atomTypes, neededParams, moleculeATs, moleculeNP,
						isDepad), ReferenceData.GEOM_WEIGHT));
	}

	private static SimpleMatrix findMockHessian(SimpleMatrix newGradient, double[] oldHessian, double[] oldGradient,
												double[] oldDir, int size) {
		SimpleMatrix s = new SimpleMatrix(oldDir);
		SimpleMatrix hessian = new SimpleMatrix(size, size);

		int count = 1;
		int index = 0;
		while (index < ((size + 1) * size) / 2) {
			for (int i = count - 1; i < size; i++) {
				hessian.set(count - 1, i, oldHessian[index]);
				hessian.set(i, count - 1, oldHessian[index]);
				index++;
			}
			count++;
		}

		SimpleMatrix y = newGradient.minus(new SimpleMatrix(oldGradient));
		double b = y.transpose().mult(s).get(0);
		SimpleMatrix A = y.mult(y.transpose()).scale(1 / b);
		double a = s.transpose().mult(hessian).mult(s).get(0);
		SimpleMatrix C = hessian.mult(s).mult(s.transpose()).mult(hessian.transpose()).scale(1 / a);

		return hessian.plus(A).minus(C);
	}
}