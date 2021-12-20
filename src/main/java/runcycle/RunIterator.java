package runcycle;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.geometry.GeometryOptimization;
import nddo.param.ParamErrorFunction;
import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import nddo.param.ParamSecondDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.optimize.ParamOptimizer;
import runcycle.optimize.ReferenceData;
import runcycle.structs.*;

import java.util.*;
import java.util.concurrent.*;

import static runcycle.State.getConverter;

public final class RunIterator implements Iterator<RunOutput>, Iterable<RunOutput> {
	private static final Logger logger = LogManager.getLogger("RunIterator.class");
	public final RunInput initialRunInput;
	private final int limit;
	private int runNumber = 0;
	private RunInput currentRunInput;
	private IMoleculeResult[] ranMolecules;

	public RunIterator(RunInput runInput) {
		this.initialRunInput = runInput;
		this.currentRunInput = runInput;
		this.limit = 0;
	}

	public RunIterator(RunInput runInput, int limit) {
		this.initialRunInput = runInput;
		this.currentRunInput = runInput;
		this.limit = limit;
	}

	public RunIterator(RunInput runInput, IMoleculeResult[] ranMolecules) {
		this.initialRunInput = runInput;
		this.currentRunInput = runInput;
		this.ranMolecules = ranMolecules;
		this.limit = 0;
	}

	public RunIterator(RunInput runInput, int limit, IMoleculeResult[] ranMolecules) {
		this.initialRunInput = runInput;
		this.currentRunInput = runInput;
		this.ranMolecules = ranMolecules;
		this.limit = limit;
	}

	public static RunOutput runOnce(RunInput ri) {
		RunIterator it = new RunIterator(ri);
		return it.next();
	}

	public static RunOutput runOnce(RunInput ri, IMoleculeResult[] ranMolecules) {
		RunIterator it = new RunIterator(ri, ranMolecules);
		return it.next();
	}

	@Override
	public boolean hasNext() {
		if (limit != 0) return runNumber < limit;
		return true;
	}

	@Override
	public RunOutput next() {
		logger.info("Run number: {}, input hash: {}", runNumber, currentRunInput.hash);

		PRun pRun = new PRun(currentRunInput);

		if (ranMolecules != null) {
			pRun.ranMolecules = ranMolecules;
			ranMolecules = null;
		}

		RunOutput output = pRun.run();

		logger.info("Run {} time taken: {}, output hash: {}\n", runNumber, output.timeTaken, output.hash);

		currentRunInput = output.getNextInput();

		runNumber++;

		return output;
	}

	@Override
	public Iterator<RunOutput> iterator() {
		return this;
	}

	public int getRunNumber() {
		return runNumber;
	}

	private static final class PRun { // stands for ParameterizationRun
		private static final Logger logger = LogManager.getLogger("PRun.class");
		private final RunInput ri;
		private final ScheduledExecutorService progressBar = Executors.newScheduledThreadPool(1);
		private IMoleculeResult[] ranMolecules;

		PRun(RunInput runInput) {
			this.ri = runInput;
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
				if (info.npMap[i] != null) newNpMap[i] = info.npMap[i].copy();
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

			return new RunOutput(ri, nextRunInfo, results.toArray(new IMoleculeResult[0]), lsw.getTime());
		}
	}

	private static final class MoleculeRun implements IMoleculeResult {
		private final NDDOAtom[] nddoAtoms, expGeom;
		private final boolean withHessian, isExpAvail;
		private final RunnableMolecule rm;
		private final double[] datum;
		private Solution s, sExp;
		private ParamGradient g;
		private ParamHessian h;
		private Atom[] newAtoms;
		private long time;

		public MoleculeRun(RunnableMolecule rm, NDDOAtom[] nddoAtoms, NDDOAtom[] expGeom, double[] datum,
						   boolean withHessian) {
			this.rm = rm;
			this.nddoAtoms = nddoAtoms;
			this.expGeom = expGeom;
			this.datum = datum;
			this.withHessian = withHessian;

			isExpAvail = expGeom != null;
		}

		public void run() {
			rm.getLogger().info("Started");
			StopWatch sw = new StopWatch();
			sw.start();


			s = GeometryOptimization.of(Solution.of(rm, nddoAtoms)).compute().getS();
			rm.getLogger().debug("Finished geometry optimization");

			if (isExpAvail) {
				sExp = Solution.of(rm, expGeom);
			}

			g = ParamGradient.of(s, datum, sExp).compute();
			rm.getLogger().debug("Finished param gradient");
			if (withHessian) h = ParamHessian.from(g).compute();
			rm.getLogger().debug("Finished param hessian");

			// stores new optimized geometry
			newAtoms = new Atom[s.atoms.length];
			for (int i = 0; i < newAtoms.length; i++) {
				newAtoms[i] = new Atom(s.atoms[i].getAtomProperties().getZ(), s.atoms[i].getCoordinates());
			}

			sw.stop();
			time = sw.getTime();


			for (int i = 1; i < 7; i++) {
				for (int j = 1; j < 7; j++) {
					System.err.println("C" + i + " C" + j);
					if (!ParamSecondDerivative.verifyEquations((SolutionR) s, 6, i, 6, j)) {
						System.err.println("clown");
						System.exit(0);
					}

					System.err.println("C" + i + " N" + j);
					if (!ParamSecondDerivative.verifyEquations((SolutionR) s, 6, i, 7, j)) {
						System.err.println("clown");
						System.exit(0);
					}

					System.err.println("N" + i + " C" + j);
					if (!ParamSecondDerivative.verifyEquations((SolutionR) s, 7, i, 6, j)) {
						System.err.println("clown");
						System.exit(0);
					}

					System.err.println("N" + i + " N" + j);
					if (!ParamSecondDerivative.verifyEquations((SolutionR) s, 7, i, 7, j)) {
						System.err.println("clown");
						System.exit(0);
					}
				}
			}
		}

		public boolean isExpAvail() {
			return isExpAvail;
		}

		public RunnableMolecule getUpdatedRm() {
			return new RunnableMolecule(rm, newAtoms, rm.expGeom, rm.datum);
		}

		public long getTime() {
			return time;
		}

		@Override
		public double getHF() {
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
		public double getGeomGradient() {
			return getE().getGeomGradient();
		}

		@Override
		public double getTotalError() {
			return getE().getTotalError();
		}

		@Override
		public double[][] getHFDerivs() {
			return getG().getHFDerivs();
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

		public ParamGradient getG() {
			return g;
		}

		public ParamHessian getH() {
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
