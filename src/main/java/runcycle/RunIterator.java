package runcycle;

import com.sun.management.OperatingSystemMXBean;
import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.geometry.GeometryOptimization;
import nddo.param.*;
import nddo.solution.Solution;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.ejml.simple.SimpleMatrix;
import runcycle.optimize.ParamOptimizer;
import runcycle.optimize.ReferenceData;
import runcycle.structs.*;

import java.lang.management.ManagementFactory;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import static runcycle.State.getConverter;

public final class RunIterator implements Iterator<RunOutput>, Iterable<RunOutput> {
	private static final OperatingSystemMXBean bean =
			(OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();

	public final RunInput initialRunInput;
	public final int limit;
	private final Logger logger;
	private final Level defaultLevel;
	private int runNumber = 0;
	private RunInput currentRunInput;
	private IMoleculeResult[] ranMolecules;

	public RunIterator(RunInput runInput) {
		this(runInput, 0, null);
	}

	public RunIterator(RunInput runInput, int limit) {
		this(runInput, limit, null);
	}

	public RunIterator(RunInput runInput, IMoleculeResult[] ranMolecules) {
		this(runInput, 0, ranMolecules);
	}

	public RunIterator(RunInput runInput, int limit, IMoleculeResult[] ranMolecules) {
		this.initialRunInput = runInput;
		this.currentRunInput = runInput;
		this.ranMolecules = ranMolecules;
		this.limit = limit;

		logger = LogManager.getLogger(initialRunInput.hash);
		defaultLevel = LogManager.getRootLogger().getLevel();
	}

	public static Atom[] getAtoms(Solution s) {
		Atom[] newAtoms = new Atom[s.atoms.length];
		for (int i = 0; i < newAtoms.length; i++) {
			newAtoms[i] = new Atom(s.atoms[i].getAtomProperties().getZ(), s.atoms[i].getCoordinates());
		}
		return newAtoms;
	}

	public RunInput getCurrentRunInput() {
		return currentRunInput;
	}

	@Override
	public boolean hasNext() {
		if (limit != 0) return runNumber < limit;
		return true;
	}

	@Override
	public RunOutput next() {
		logger.info("Run number: {}, input hash: {}", runNumber, currentRunInput.hash);

		PRun pRun = new PRun(currentRunInput, runNumber, logger);

		if (ranMolecules != null) {
			pRun.ranMolecules = ranMolecules;
			ranMolecules = null;
		}

		Configurator.setRootLevel(defaultLevel);

		try {
			RunOutput output = pRun.run();

			logger.info("Run {} time taken: {}, output hash: {}", runNumber, output.timeTaken, output.hash);

			currentRunInput = output.nextInput;

			runNumber++;

			return output;
		} catch (Exception e) {
			logger.error("{} errored!", runNumber, e);
			throw e;
		}
	}

	@Override
	public Iterator<RunOutput> iterator() {
		return this;
	}

	public int getRunNumber() {
		return runNumber;
	}

	private static final class PRun { // stands for ParameterizationRun
		private final Logger logger;
		private final RunInput ri;
		private final int runNumber;
		private final ScheduledExecutorService progressBar = Executors.newScheduledThreadPool(1);
		private IMoleculeResult[] ranMolecules;

		PRun(RunInput runInput, int runNumber, Logger logger) {
			this.ri = runInput;
			this.runNumber = runNumber;
			this.logger = logger;
		}

		private static int getMaxMoleculeIndex(RunnableMolecule[] rms) {
			int max = 0;

			for (RunnableMolecule rm : rms) {
				if (rm.index > max) max = rm.index;
			}

			return max;
		}

		public RunOutput run() {
			final InputInfo info = ri.info;
			final RunnableMolecule[] rms = ri.molecules;

			boolean[] isDones = new boolean[getMaxMoleculeIndex(rms) + 1];

			StopWatch lsw = new StopWatch();
			lsw.start();

			if (logger.isInfoEnabled()) {
				AtomicInteger lastLeftCount = new AtomicInteger(ri.molecules.length);
				AtomicInteger count = new AtomicInteger(0);
				Runnable mLeft = () -> {
					List<String> left = new ArrayList<>();
					int doneCount = 0;

					for (RunnableMolecule rm : rms) {
						if (!isDones[rm.index]) left.add(rm.debugName());
						else doneCount++;
					}

					int leftCount = left.size();

					if (count.get() % 6 == 1) {
						if (lastLeftCount.get() - leftCount == 0) {
							logger.warn("Stubborn molecule(s) detected, increasing log level...");
							Configurator.setRootLevel(Level.TRACE);
						}
						lastLeftCount.set(leftCount);
					}

					count.incrementAndGet();

					Collections.sort(left);

					int totalCount = doneCount + leftCount;
					String percent = String.format("%.2f", 100.0 * doneCount / totalCount);

					long time = lsw.getTime(TimeUnit.SECONDS);
					double systemCpuLoad = bean.getSystemCpuLoad();
					double eta = time * (1.0 * totalCount / doneCount);
					logger.info("Run {}, time: {} s, CPU load: {}, ETA: {} s, {}/{} left ({}% done): {}",
							runNumber, time, String.format("%.2f", systemCpuLoad), String.format("%.2f", eta),
							leftCount, totalCount, percent, left);
				};

				int wait = 30;
				progressBar.scheduleAtFixedRate(mLeft, wait, wait, TimeUnit.SECONDS);
			}


			// create tasks to run in parallel and then runs them
			// excludes ran molecules from previous dynamic output
			List<RecursiveTask<MoleculeRun>> moleculeTasks = new ArrayList<>();
			for (RunnableMolecule rm : rms) {
				boolean isDoneAlready = false;

				if (ranMolecules != null) {
					for (IMoleculeResult mo : ranMolecules) {
						if (mo.getUpdatedRm().getIndex() == rm.index) {
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
			Collections.shuffle(moleculeTasks, new Random(Constants.RANDOM_SEED));
			for (ForkJoinTask<MoleculeRun> task : ForkJoinTask.invokeAll(moleculeTasks)) {
				results.add(task.join());
			}

			results.sort(Comparator.comparingInt(x -> x.getUpdatedRm().getIndex()));


			// combined length of all differentiated params
			int paramLength = 0;
			for (int[] param : info.neededParams) paramLength += param.length;

			// processing results
			ParamOptimizer opt = new ParamOptimizer();
			double ttError = 0;
			double[] ttGradient = new double[paramLength];
			double[][] ttHessian = new double[paramLength][paramLength];

			for (IMoleculeResult result : results) {
				int[] moleculeATs = result.getUpdatedRm().getMats();
				int[][] moleculeNPs = result.getUpdatedRm().getMnps();
				boolean isDepad = true;

				ttError += result.getTotalError();

				// add things to ParamOptimizer
				double[] datum = result.getUpdatedRm().getDatum();

				opt.addData(new ReferenceData(datum[0], result.getHf(),
						ParamGradient.combine(result.getHfDerivs(), info.atomTypes, info.neededParams,
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
					opt.addData(new ReferenceData(0, result.getGeomGradMag(),
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


				double[][] h = ParamHessian.padHessian(result.getHessian(), result.getUpdatedRm().getMats(),
						info.atomTypes, info.neededParams);

				boolean hasNan = false;
				for (int i = 0; i < h.length; i++) {
					for (int j = 0; j < h[0].length; j++) {
						if (Double.isNaN(h[i][j])) {
							hasNan = true;
						}
						else ttHessian[i][j] += h[i][j];
					}
				}

				if (hasNan) {
					logger.warn("NaN in Hessian! {}: \n{}\n{}", result.getUpdatedRm().debugName(), g,
							Arrays.deepToString(h));
				}
			}

			// optimizes params based on this run and gets new search direction
			SimpleMatrix newGradient = new SimpleMatrix(ttGradient);
			SimpleMatrix newHessian = new SimpleMatrix(ttHessian);

			double[] dir = opt.optimize(newHessian, newGradient);

			// generating nextRunInfo
			NDDOParams[] newNpMap = new NDDOParams[info.npMap.length];
			for (int i = 0; i < newNpMap.length; i++) {
				if (info.npMap[i] != null) newNpMap[i] = info.npMap[i].copy();
			}

			int n = 0;
			for (int atomI = 0; atomI < info.atomTypes.length; atomI++) {
				for (int neededParam : info.neededParams[atomI]) {
					newNpMap[info.atomTypes[atomI]].modifyParam(neededParam, dir[n]);
					n++;
				}
			}


			IMoleculeResult[] resultsArray = results.toArray(new IMoleculeResult[0]);
			InputInfo nextRunInfo = new InputInfo(info.atomTypes, info.neededParams, newNpMap);
			RunnableMolecule[] nextRunRms = new RunnableMolecule[resultsArray.length];

			for (int i = 0; i < nextRunRms.length; i++) {
				nextRunRms[i] = (RunnableMolecule) resultsArray[i].getUpdatedRm();
			}

			RunInput nextInput = new RunInput(nextRunInfo, nextRunRms);


			lsw.stop();
			progressBar.shutdownNow();

			logger.info("Total error: {}", ttError);

			return new RunOutput(resultsArray, lsw.getTime(), ttError, ttGradient, ttHessian, ri, nextInput);
		}
	}

	public static final class MoleculeRun implements IMoleculeResult {
		private boolean withHessian, isExpAvail;
		private RunnableMolecule rm, updatedRm;
		private Solution s;
		private IParamGradient g;
		private IParamHessian h;
		private long time;

		public MoleculeRun(RunnableMolecule rm, NDDOAtom[] nddoAtoms, NDDOAtom[] expGeom, double[] datum,
						   boolean withHessian) {
			try {
				this.rm = rm;
				this.withHessian = withHessian;
				this.isExpAvail = expGeom != null;

				rm.getLogger().debug("Started");
				StopWatch sw = StopWatch.createStarted();

				Level orig = rm.getLogger().getLevel();
				ScheduledExecutorService logIncreaser = Executors.newScheduledThreadPool(1);
				logIncreaser.schedule(() -> {
					if (sw.getTime() >= nddo.State.config.log_increase_time) {
						rm.getLogger().warn("Stubborn molecule detected- increasing log level...");
						Configurator.setLevel(rm.getLogger().getName(), Level.TRACE);
					}
				}, nddo.State.config.log_increase_time, TimeUnit.SECONDS);

				Solution initialS = Solution.of(rm, nddoAtoms, rm.densityMatrices);
				rm.getLogger().debug("Finished initial solution computation");

				s = GeometryOptimization.of(initialS).compute().getS();
				rm.getLogger().debug("Finished geometry optimization");

				Solution sExp = isExpAvail ? Solution.of(rm, expGeom, rm.densityMatricesExp) : null;
				g = new ParamGradientNew(s, datum, sExp);
				rm.getLogger().debug("Finished param gradient");

				h = withHessian ? new ParamHessianNew((ParamGradientNew) g) : null;
				rm.getLogger().debug("Finished param hessian");

				rm.densityMatrices = rm.restricted ?
						new double[][]{s.densityMatrix().getDDRM().data} :
						new double[][]{s.alphaDensity().getDDRM().data, s.betaDensity().getDDRM().data};

				if (isExpAvail) rm.densityMatricesExp = rm.restricted ?
						new double[][]{sExp.densityMatrix().getDDRM().data} :
						new double[][]{sExp.alphaDensity().getDDRM().data, sExp.betaDensity().getDDRM().data};

				// stores new optimized geometry
				updatedRm = new RunnableMolecule(rm, getAtoms(s), rm.expGeom, rm.datum);

				sw.stop();
				time = sw.getTime();

				logIncreaser.shutdownNow();
				Configurator.setLevel(rm.getLogger().getName(), orig);

				rm.getLogger().info("Finished in {}", String.format("%06d", time));
			} catch (Exception e) {
				e.printStackTrace();
				rm.getLogger().error(e);
				System.exit(1);
			}
		}

		public boolean isExpAvail() {
			return isExpAvail;
		}

		public RunnableMolecule getUpdatedRm() {
			return updatedRm;
		}

		public long getTime() {
			return time;
		}

		@Override
		public double getHf() {
			return getS().hf;
		}

		@Override
		public double getDipole() {
			return getS().dipole;
		}

		@Override
		public double getIE() {
			return -getS().homo;
		}

		@Override
		public double getGeomGradMag() {
			return getE().getGeomGradMag();
		}

		@Override
		public double getTotalError() {
			return getE().getTotalError();
		}

		@Override
		public double[][] getHfDerivs() {
			return getG().getHfDerivs();
		}

		@Override
		public double[][] getDipoleDerivs() {
			return getG().getDipoleDerivs();
		}

		@Override
		public double[][] getIEDerivs() {
			return getG().getIEDerivs();
		}

		@Override
		public double[][] getGeomDerivs() {
			return getG().getGeomDerivs();
		}

		@Override
		public double[][] getTotalGradients() {
			return getG().getTotalGradients();
		}

		@Override
		public double[][] getHessian() {
			if (withHessian) return h.getHessian();
			else throw new IllegalStateException("Hessian not found for molecule: " + rm.debugName());
		}

		public IParamGradient getG() {
			return g;
		}

		public IParamHessian getH() {
			return h;
		}

		/**
		 * @return Original, unoptimized Solution object.
		 */
		public Solution getS() {
			return s;
		}

		public ParamErrorFunction getE() {
			return g.getE();
		}
	}
}
