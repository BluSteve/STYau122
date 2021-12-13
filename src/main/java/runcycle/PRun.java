package runcycle;

import nddo.NDDOAtom;
import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.input.InputInfo;
import runcycle.input.RawParams;
import runcycle.optimize.ParamOptimizer;
import runcycle.optimize.ReferenceData;
import runcycle.output.MoleculeOutput;
import runcycle.output.RunOutput;
import runcycle.structs.RunnableMolecule;

import java.util.*;
import java.util.concurrent.*;

import static runcycle.State.getConverter;

public class PRun { // stands for ParameterizationRun
	private static final Logger logger = LogManager.getLogger();
	public final InputInfo info;
	public final RunnableMolecule[] rms;
	private final ScheduledExecutorService progressBar = Executors.newScheduledThreadPool(1);
	private MoleculeOutput[] ranMolecules;

	public PRun(InputInfo info, RunnableMolecule[] rms) {
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

	public void setRanMolecules(MoleculeOutput[] ranMolecules) {
		this.ranMolecules = ranMolecules;
	}

	public RunOutput run(boolean isRunHessian) {
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
				for (MoleculeOutput mo : ranMolecules) {
					if (mo.runnableMolecule.index == rm.index) {
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
						NDDOAtom[] atoms = getConverter().convert(rm.atoms);
						NDDOAtom[] expGeom = rm.expGeom == null ? null : getConverter().convert(rm.expGeom);

						MoleculeRun mr = new MoleculeRun(rm, atoms, expGeom, isRunHessian);
						mr.run();

						isDones[rm.index] = true;
						return mr;
					}
				});
			}
		}


		// adds previously ran molecules into the final results list
		List<MoleculeResult> results = new ArrayList<>(rms.length);

		if (ranMolecules != null) {
			for (MoleculeOutput mr : ranMolecules) {
				results.add(new MoleculeRan(mr));
			}
		}

		// shuffles run requests and runs them, then deshuffles them
		Collections.shuffle(moleculeTasks, new Random(123));
		for (ForkJoinTask<MoleculeRun> task : ForkJoinTask.invokeAll(moleculeTasks)) {
			results.add(task.join());
		}
		results.sort(Comparator.comparingInt(x -> x.getRm().index));


		// processing results
		ParamOptimizer opt = new ParamOptimizer();
		double[] ttGradient = new double[paramLength];
		double[][] ttHessian = isRunHessian ? new double[paramLength][paramLength] : null;
		MoleculeOutput[] mos = new MoleculeOutput[results.size()];

		int index = 0;
		for (MoleculeResult result : results) {
			int[] moleculeATs = result.getRm().mats;
			int[][] moleculeNPs = result.getRm().mnps;
			boolean isDepad = true;


			// add things to ParamOptimizer
			double[] datum = result.getRm().datum;

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


			if (isRunHessian) {
				double[][] h = ParamHessian.padHessian(result.getHessian(), result.getRm().mats, info.atomTypes,
						info.neededParams);

				for (int i = 0; i < h.length; i++) {
					for (int j = 0; j < h[0].length; j++)
						ttHessian[i][j] += h[i][j];
				}
			}


			mos[index] = MoleculeOutput.from(result);
			index++;
		}

		// optimizes params based on this run and gets new search direction
		SimpleMatrix newGradient = new SimpleMatrix(ttGradient);
		SimpleMatrix newHessian = isRunHessian ? new SimpleMatrix(ttHessian) :
				findMockHessian(newGradient, info.params.lastHessian, info.params.lastGradient, info.params.lastDir,
						paramLength);

		double[] dir = opt.optimize(newHessian, newGradient);

		// generating nextRunInfo
		double[][] nddoParams = new double[info.params.nddoParams.length][]; // deepcopy nddoParams

		for (int i = 0; i < info.params.nddoParams.length; i++) {
			nddoParams[i] = info.params.nddoParams[i].clone();
		}

		int n = 0;
		for (int atomI = 0; atomI < info.neededParams.length; atomI++) {
			for (int paramI : info.neededParams[atomI]) {
				nddoParams[atomI][paramI] += dir[n];
				n++;
			}
		}

		RawParams rawParams = new RawParams(nddoParams, dir, ttGradient, ParamHessian.utify(newHessian.toArray2()));
		InputInfo nextRunInfo = new InputInfo(info.atomTypes, info.neededParams, rawParams);

		return new ConcreteRunOutput(nextRunInfo, mos);
	}

	private static class ConcreteRunOutput implements RunOutput {
		InputInfo nextRunInfo;
		MoleculeOutput[] mos;

		public ConcreteRunOutput(InputInfo nextRunInfo, MoleculeOutput[] mos) {
			this.nextRunInfo = nextRunInfo;
			this.mos = mos;
		}

		@Override
		public InputInfo getNextRunInfo() {
			return nextRunInfo;
		}

		@Override
		public MoleculeOutput[] getMOs() {
			return mos;
		}
	}
}
