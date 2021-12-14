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
import runcycle.structs.Params;
import runcycle.structs.RunOutput;
import runcycle.structs.RunnableMolecule;

import java.util.*;
import java.util.concurrent.*;

import static runcycle.State.getConverter;

public class RunIterator implements Iterator<RunOutput> {
	private static final Logger logger = LogManager.getLogger();
	private final InputInfo initialInfo;
	private final RunnableMolecule[] initialRms;
	private int runNumber = 0, limit = 0;
	private InputInfo currentInfo;
	private RunnableMolecule[] currentRms;
	private boolean withHessian;

	public RunIterator(InputInfo initialInfo, RunnableMolecule[] initialRms, boolean withHessian) {
		this.initialInfo = initialInfo;
		this.initialRms = initialRms;
		this.withHessian = withHessian;

		this.currentInfo = initialInfo;
		this.currentRms = initialRms;
	}

	@Override
	public boolean hasNext() {
		if (limit != 0) return runNumber < limit;
		return true;
	}

	@Override
	public RunOutput next() {
		logger.info("Run number: {}, hessian: {}", runNumber, withHessian);

		RunOutput output = new PRun(currentInfo, currentRms).run(withHessian);
		runNumber++;

		currentInfo = output.nextRunInfo;

		currentRms = new RunnableMolecule[output.results.length];
		for (int i = 0; i < currentRms.length; i++) {
			currentRms[i] = output.results[i].getUpdatedRm();
		}

		logger.info("Run {} time taken: {}", runNumber, output.timeTaken);

		return output;
	}

	public InputInfo getCurrentInfo() {
		return currentInfo;
	}

	public RunnableMolecule[] getCurrentRms() {
		return currentRms;
	}

	public boolean isWithHessian() {
		return withHessian;
	}

	public void setWithHessian(boolean withHessian) {
		this.withHessian = withHessian;
	}

	public int getLimit() {
		return limit;
	}

	public void setLimit(int limit) {
		this.limit = limit;
	}

	public InputInfo getInitialInfo() {
		return initialInfo;
	}

	public RunnableMolecule[] getInitialRms() {
		return initialRms;
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

		PRun(InputInfo info, RunnableMolecule[] rms) {
			this.info = info;
			this.rms = rms;
		}

		private static int getMaxMoleculeIndex(RunnableMolecule[] rms) {
			int max = 0;

			for (RunnableMolecule rm : rms) {
				if (rm.index > max) max = rm.index;
			}

			return max;
		}

		private static SimpleMatrix findMockHessian(SimpleMatrix newGradient, double[] oldHessian,
													double[] oldGradient,
													double[] oldDir, int size) {
			SimpleMatrix s = new SimpleMatrix(oldDir);
			SimpleMatrix hessian = new SimpleMatrix(size, size);

			int count = 1;
			int index = 0;
			while (index < (size + 1) * size / 2) {
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

		public InputInfo getInfo() {
			return info;
		}

		public RunnableMolecule[] getRms() {
			return rms;
		}

		public RunOutput run(boolean withHessian) {
			boolean[] isDones = new boolean[getMaxMoleculeIndex(rms)];

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
							NDDOAtom[] atoms = getConverter().convert(rm.atoms, info.params.npMap);
							NDDOAtom[] expGeom = rm.expGeom == null ?
									null : getConverter().convert(rm.expGeom, info.params.npMap);

							MoleculeRun mr = new MoleculeRun(rm, atoms, expGeom, rm.datum, withHessian);
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
			double[] ttGradient = new double[paramLength];
			double[][] ttHessian = withHessian ? new double[paramLength][paramLength] : null;

			for (IMoleculeResult result : results) {
				int[] moleculeATs = result.getUpdatedRm().mats;
				int[][] moleculeNPs = result.getUpdatedRm().mnps;
				boolean isDepad = true;


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


				if (withHessian) {
					double[][] h =
							ParamHessian.padHessian(result.getHessian(), result.getUpdatedRm().mats, info.atomTypes,
									info.neededParams);

					for (int i = 0; i < h.length; i++) {
						for (int j = 0; j < h[0].length; j++)
							ttHessian[i][j] += h[i][j];
					}
				}
			}

			// optimizes params based on this run and gets new search direction
			SimpleMatrix newGradient = new SimpleMatrix(ttGradient);
			SimpleMatrix newHessian = withHessian ? new SimpleMatrix(ttHessian) :
					findMockHessian(newGradient, info.params.lastHessian, info.params.lastGradient,
							info.params.lastDir,
							paramLength);

			double[] dir = opt.optimize(newHessian, newGradient);

			// generating nextRunInfo
			NDDOParams[] newNpMap = new NDDOParams[info.params.npMap.length];
			for (int i = 0; i < newNpMap.length; i++) {
				newNpMap[i] = info.params.npMap[i].clone();
			}

			int n = 0;
			for (int atomI = 0; atomI < info.atomTypes.length; atomI++) {
				for (int neededParam : info.neededParams[atomI]) {
					newNpMap[info.atomTypes[atomI]].modifyParam(neededParam, dir[n]);
					n++;
				}
			}

			Params params = new Params(newNpMap, dir, ttGradient, ParamHessian.utify(newHessian.toArray2()));
			InputInfo nextRunInfo = new InputInfo(info.atomTypes, info.neededParams, params);

			lsw.stop();
			progressBar.shutdownNow();

			return new RunOutput(nextRunInfo, results.toArray(new IMoleculeResult[0]), lsw.getTime());
		}
	}
}
