package runcycle;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.optimize.ParamOptimizer;
import runcycle.optimize.ReferenceData;
import runcycle.structs.InputInfo;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;
import runcycle.structs.RunnableMolecule;

import java.util.*;
import java.util.concurrent.*;

import static runcycle.State.getConverter;

public class RunIterator implements Iterator<RunOutput> {
	private static final Logger logger = LogManager.getLogger();
	private final RunInput initialRunInput;
	private int runNumber = 0, limit = 0;
	private RunInput currentRunInput;

	public RunIterator(RunInput runInput) {
		this.initialRunInput = runInput;

		this.currentRunInput = runInput;
	}

	@Override
	public boolean hasNext() {
		if (limit != 0) return runNumber < limit;
		return true;
	}

	@Override
	public RunOutput next() {
		logger.info("Run number: {}", runNumber);

		RunOutput output = new PRun(currentRunInput).run();

		RunnableMolecule[] currentRms = new RunnableMolecule[output.results.length];
		for (int i = 0; i < currentRms.length; i++) {
			currentRms[i] = output.results[i].getUpdatedRm();
		}

		currentRunInput = new RunInput(output.nextRunInfo, currentRms);

		logger.info("Run {} time taken: {}\n", runNumber, output.timeTaken);

		runNumber++;

		return output;
	}

	public RunInput getInitialRunInput() {
		return initialRunInput;
	}

	public RunInput getCurrentRunInput() {
		return currentRunInput;
	}

	public int getLimit() {
		return limit;
	}

	public void setLimit(int limit) {
		this.limit = limit;
	}

	public int getRunNumber() {
		return runNumber;
	}

	private static class PRun { // stands for ParameterizationRun
		private static final Logger logger = LogManager.getLogger();
		private final InputInfo info;
		private final RunnableMolecule[] rms;
		private final ScheduledExecutorService progressBar = Executors.newScheduledThreadPool(1);
		private IMoleculeResult[] ranMolecules;

		PRun(RunInput runInput) {
			this.info = runInput.info;
			this.rms = runInput.molecules;
		}

		private static int getMaxMoleculeIndex(RunnableMolecule[] rms) {
			int max = 0;

			for (RunnableMolecule rm : rms) {
				if (rm.index > max) max = rm.index;
			}

			return max;
		}

		public InputInfo getInfo() {
			return info;
		}

		public RunnableMolecule[] getRms() {
			return rms;
		}

		public RunOutput run() {
			boolean[] isDones = new boolean[getMaxMoleculeIndex(rms) + 1];

			if (logger.isInfoEnabled()) {
				Runnable mLeft = () -> {
					List<String> done = new ArrayList<>();
					List<String> left = new ArrayList<>();

					for (RunnableMolecule rm : rms) {
						if (isDones[rm.index]) done.add(rm.debugName());
						else left.add(rm.debugName());
					}

					done.sort(String::compareTo);
					left.sort(String::compareTo);

					logger.info("Molecules done ({}): {}", done.size(), done);
					logger.info("Molecules left ({}): {}", left.size(), left);
				};

				int wait = 10;
				progressBar.scheduleAtFixedRate(mLeft, wait, wait, TimeUnit.SECONDS);
			}

			StopWatch lsw = new StopWatch();
			lsw.start();


			// combined length of all differentiated params
			int paramLength = 0;
			for (int[] param : info.neededParams) paramLength += param.length;


			// create tasks to run in parallel and then runs them
			// excludes ran molecules from previous dynamic output
			List<RecursiveTask<MoleculeRun>> moleculeTasks = new ArrayList<>();
			for (RunnableMolecule rm : rms) {
				boolean isDoneAlready = false;

				if (ranMolecules != null) {
					for (IMoleculeResult mo : ranMolecules) {
						if (mo.getUpdatedRm().index == rm.index) {
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
							NDDOAtom[] atoms = getConverter().convert(rm.atoms, info.npMap);
							NDDOAtom[] expGeom = rm.expGeom == null ?
									null : getConverter().convert(rm.expGeom, info.npMap);

							MoleculeRun mr = new MoleculeRun(rm, atoms, expGeom, rm.datum, true);
							mr.run();

							isDones[rm.index] = true;
							return mr;
						}
					});
				}
			}


			// adds previously ran molecules into the final results list
			List<IMoleculeResult> results = new ArrayList<>(rms.length);

			if (ranMolecules != null) results.addAll(List.of(ranMolecules));

			// shuffles run requests and runs them, then deshuffles them
			Collections.shuffle(moleculeTasks, new Random(123));
			for (ForkJoinTask<MoleculeRun> task : ForkJoinTask.invokeAll(moleculeTasks)) {
				results.add(task.join());
			}
			results.sort(Comparator.comparingInt(x -> x.getUpdatedRm().index));


			// processing results
			ParamOptimizer opt = new ParamOptimizer();
			double ttError = 0;
			double[] ttGradient = new double[paramLength];
			double[][] ttHessian = new double[paramLength][paramLength];

			for (IMoleculeResult result : results) {
				int[] moleculeATs = result.getUpdatedRm().mats;
				int[][] moleculeNPs = result.getUpdatedRm().mnps;
				boolean isDepad = true;

				ttError += result.getTotalError();

				// add things to ParamOptimizer
				double[] datum = result.getUpdatedRm().datum;

				opt.addData(new ReferenceData(datum[0], result.getHF(),
						ParamGradient.combine(result.getHFDerivs(), info.atomTypes, info.neededParams,
								moleculeATs, moleculeNPs, isDepad),
						ReferenceData.HF_WEIGHT));

				if (datum[1] != 0) {
					opt.addData(new ReferenceData(datum[1], result.getDipole(),
							ParamGradient.combine(result.getDipoleDerivs(), info.atomTypes, info.neededParams,
									moleculeATs, moleculeNPs, isDepad),
							ReferenceData.DIPOLE_WEIGHT));
				}

				if (datum[2] != 0) {
					opt.addData(new ReferenceData(datum[2], result.getIE(),
							ParamGradient.combine(result.getIEDerivs(), info.atomTypes, info.neededParams,
									moleculeATs, moleculeNPs, isDepad),
							ReferenceData.IE_WEIGHT));
				}

				if (result.isExpAvail()) {
					opt.addData(new ReferenceData(0, result.getGeomGradient(),
							ParamGradient.combine(result.getGeomDerivs(), info.atomTypes, info.neededParams,
									moleculeATs, moleculeNPs, isDepad),
							ReferenceData.GEOM_WEIGHT));
				}


				// ttGradient is sum of totalGradients across molecules
				double[] g = ParamGradient.combine(result.getTotalGradients(), info.atomTypes, info.neededParams,
						moleculeATs, moleculeNPs, isDepad);

				for (int i = 0; i < g.length; i++) {
					ttGradient[i] += g[i];
				}


				double[][] h = ParamHessian.padHessian(result.getHessian(), result.getUpdatedRm().mats,
						info.atomTypes, info.neededParams);

				for (int i = 0; i < h.length; i++) {
					for (int j = 0; j < h[0].length; j++)
						ttHessian[i][j] += h[i][j];
				}
			}

			// optimizes params based on this run and gets new search direction
			SimpleMatrix newGradient = new SimpleMatrix(ttGradient);
			SimpleMatrix newHessian = new SimpleMatrix(ttHessian);

			double[] dir = opt.optimize(newHessian, newGradient);

			// generating nextRunInfo
			NDDOParams[] newNpMap = new NDDOParams[info.npMap.length];
			for (int i = 0; i < newNpMap.length; i++) {
				if (info.npMap[i] != null) newNpMap[i] = info.npMap[i].clone();
			}

			int n = 0;
			for (int atomI = 0; atomI < info.atomTypes.length; atomI++) {
				for (int neededParam : info.neededParams[atomI]) {
					newNpMap[info.atomTypes[atomI]].modifyParam(neededParam, dir[n]);
					n++;
				}
			}

			InputInfo nextRunInfo = new InputInfo(info.atomTypes, info.neededParams, newNpMap);

			lsw.stop();
			progressBar.shutdownNow();

			logger.info("Total error: {}", ttError);

			return new RunOutput(nextRunInfo, results.toArray(new IMoleculeResult[0]), lsw.getTime());
		}
	}
}
