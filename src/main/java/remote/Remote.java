package remote;

import frontend.FrontendConfig;
import frontend.JsonIO;
import host.NodeDisconnectedException;
import host.Server;
import nddo.NDDOParams;
import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.IMoleculeResult;
import runcycle.optimize.IParamOptimizer;
import runcycle.optimize.LLAOptimizer;
import runcycle.optimize.ReferenceData;
import runcycle.structs.*;
import tools.Utils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Queue;
import java.util.Scanner;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static frontend.FrontendConfig.config;

public class Remote {
	private static double[] timeTaken;

	// todo make websocket faster by using java serialization vs gson
	public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
		FrontendConfig.init();
		Logger logger = LogManager.getLogger();

		logger.info("Date compiled: {}", Utils.getResource("version.txt"));

		Files.createDirectories(Path.of("pastinputs"));
		Files.createDirectories(Path.of("outputs"));
		Files.createDirectories(Path.of("summaries"));

		Server<AdvancedMachine> server = new Server<>(56701, AdvancedMachine.class);
		server.setOnConnectListener(AdvancedMachine::init);
		while (server.getMachineCount() == 0) {
		}

		Scanner s = new Scanner(System.in);
		while (!s.nextLine().equals("ok")) {
		}

		AdvancedMachine[] machines = server.getMachines();
		Arrays.sort(machines);
		timeTaken = new double[machines.length];
		logger.info("Machines ({}): {}", machines.length, Arrays.toString(machines));


//		 build initial RunInput object
//		String pnFile = null;
//		String pFile = Files.readString(Path.of("params.csv"));
//		String mFile = Files.readString(Path.of("molecules.txt"));
//
//		AdvancedMachine bestMachine = machines[0];
//		RunInput runInput = bestMachine.buildMolecules(pnFile, pFile, mFile);
//		JsonIO.write(runInput, "remote-input");
//		logger.info("Finished initializing molecules.");

		RunInput runInput = JsonIO.readInput("remote-input");

		for (int i = 0; i < runInput.molecules.length; i++) {
			if (runInput.molecules[i].index != i) logger.warn("index mismatch: " + i);
		}

		// creating endingIndices to group molecules by
		int length = runInput.molecules.length;
		int[] endingIndices = getEndingIndices(machines, length);

		logger.info("Length: {}, ending indices: {}, machines: {}", length, endingIndices, machines);


		// do runs
		int i = config.starting_run_num;
		try {
			RunInput currentRunInput = runInput;
			while (i < config.num_runs) {
				JsonIO.writeAsync(currentRunInput, String.format("pastinputs/%04d-%s", i, currentRunInput.hash));

				if (machines.length != server.getMachineCount()) {
					machines = server.getMachines();
					endingIndices = getEndingIndices(machines, length);
					timeTaken = new double[machines.length];
				}

				logger.info("Machines ({}): {}", machines.length, Arrays.toString(machines));
				logger.info("Run number: {}, input hash: {}", i, currentRunInput.hash);

				RunOutput ro;
				try {
					ro = run(machines, endingIndices, currentRunInput);
				} catch (NodeDisconnectedException e) {
					logger.warn("{} disconnected; restarting run...", e.getMachine().name);

					server.remove(e.getMachine().name);

					if (server.getMachineCount() == 0) {
						throw new RuntimeException("No more machines left to use! Quitting...");
					}

					continue;
				}

				JsonIO.writeAsync(ro, String.format("outputs/%04d-%s-%s", i, ro.inputHash, ro.hash));

				ROSummary summary = getSummary(ro);
				JsonIO.write(summary, String.format("summaries/%04d-%s-%s", i, ro.inputHash, ro.hash));

				logger.info("Run {} time taken: {}, output hash: {}, next input hash: {}\n\n", i, ro.timeTaken,
						ro.hash, ro.nextInputHash);

				currentRunInput = ro.nextInput;

				double max = 0;
				for (double v : timeTaken) if (v > max) max = v;
				for (int j = 0; j < timeTaken.length; j++) timeTaken[j] /= max;

				double spread = Utils.sd(timeTaken);
				if (spread > config.reconf_power_threshold) {
					logger.info("Spread = {}, recalibrating power of machines... (curr={})", spread,
							endingIndices);

					for (int j = 0; j < machines.length; j++) {
						machines[j].power =
								(machines[j].power + (endingIndices[j + 1] - endingIndices[j]) / timeTaken[j]) / 2;
					}
					endingIndices = getEndingIndices(machines, length);

					logger.info("Finished recalibrating power of machines (new={})", endingIndices);

					Arrays.stream(machines).parallel().forEach(AdvancedMachine::updatePower);

					logger.info("Uploaded new powers: {}\n\n", Arrays.toString(machines));
				}
				else {
					logger.info("Spread = {}", spread);
				}

				timeTaken = new double[machines.length];

				i++;
			}
		} catch (Exception e) {
			logger.error("Run {} errored!", i, e);
			throw e;
		}

		System.exit(0);
	}

	private static ROSummary getSummary(RunOutput ro) {
		ROSummary summary = new ROSummary();

		summary.error = ro.ttError;
		summary.lambda = ro.finalLambda;
		summary.hessian = ro.ttHessian;
		summary.evs = Utils.symEigen(new SimpleMatrix(summary.hessian))[1].diag().getDDRM().data;
		summary.params = ro.nextInput.info.getParams();

		return summary;
	}

	private static class ROSummary {
		double error, lambda;
		double[] evs;
		double[][] params;
		double[][] hessian;
	}

	private static int[] getEndingIndices(AdvancedMachine[] machines, int length) {
		double totalpower = 0;
		for (AdvancedMachine executor : machines) totalpower += executor.power;

		int[] endingIndices = new int[machines.length + 1];
		for (int i = 0; i < machines.length; i++) {
			int inte = endingIndices[i]; // previous index
			int end = (int) Math.round(length * (machines[i].power / totalpower) + inte);
			endingIndices[i + 1] = end;
		}
		endingIndices[endingIndices.length - 1] = length; // just in case rounding issue

		return endingIndices;
	}

	public static RunOutput run(AdvancedMachine[] machines, int[] endingIndices, RunInput runInput) {
		Logger logger = LogManager.getLogger(runInput.hash);

		RunnableMolecule[] rms = runInput.molecules;
		InputInfo info = runInput.info;
		int nMachines = machines.length;

		StopWatch sw = StopWatch.createStarted();

		// grouping molecules based on coreCount
		Utils.shuffleArray(rms); // shuffles whole rms, should be the same across runs due to same seed
		RunnableMolecule[][] rms2d = new RunnableMolecule[nMachines][];
		for (int i = 1; i < endingIndices.length; i++) {
			RunnableMolecule[] rmsubset = Arrays.copyOfRange(rms, endingIndices[i - 1], endingIndices[i]);
			Arrays.sort(rmsubset, Comparator.comparingInt(rm -> rm.index)); // sorted
			rms2d[i - 1] = rmsubset;
		}
		Arrays.sort(rms, Comparator.comparingInt(rm -> rm.index)); // deshuffles rms


		// doing the actual computation
		IMoleculeResult[][] results2d = new IMoleculeResult[nMachines][];
		AtomicInteger doneCount = new AtomicInteger(0);
		AtomicInteger doneMachineCount = new AtomicInteger(0);
		Queue<AdvancedMachine> errored = new ConcurrentLinkedQueue<>();
		IntStream.range(0, rms2d.length).parallel().forEach(i -> { // multithreaded uploading
			AdvancedMachine machine = machines[i];
			Logger machineLogger = machine.logger;

			try {
				StopWatch sw2 = StopWatch.createStarted();
				AdvancedMachine.SubsetResult result = machine.runMolecules(rms2d[i], info, runInput.hash);
				sw2.stop();

				timeTaken[i] += sw2.getTime(); // inclusive of latency

				results2d[i] = result.results;
				Arrays.sort(results2d[i], Comparator.comparingInt(r -> r.getUpdatedRm().index));

				if (results2d[i].length != rms2d[i].length) {
					throw new RuntimeException(machine.ip + ": " + rms2d[i].length + " " + results2d[i].length);
				}

				if (machineLogger.isInfoEnabled()) machineLogger.info("\n" + machine.getLogs());

				logger.info("{} molecules finished from {}/{} machines", doneCount.addAndGet(results2d[i].length),
						doneMachineCount.incrementAndGet(), nMachines);
			} catch (Error e) {
				machineLogger.error("", e);
				errored.add(machine);
			}
		});

		if (errored.size() > 0) {
			throw new NodeDisconnectedException(errored.peek());
		}

		IMoleculeResult[] results = Stream.of(results2d).flatMap(Arrays::stream).sorted(
				Comparator.comparingInt(r -> r.getUpdatedRm().index)).toArray(IMoleculeResult[]::new);

		for (int i = 0; i < results.length; i++) {
			for (int j = i + 1; j < results.length; j++) {
				if (results[i].getUpdatedRm().index == results[j].getUpdatedRm().index) {
					throw new RuntimeException("Duplicate found! " + results[i].getUpdatedRm().index);
				}
			}
		}

		// processing results
		int paramLength = 0; // combined length of all differentiated params
		for (int[] param : info.neededParams) paramLength += param.length;


		LastRunInfo lri;
		// todo make abstraction barrier higher for paramopt. put tt hessian and gradient in PO
		LLAOptimizer lla = new LLAOptimizer();
		lri = null;

		// optimizes params based on this run
		double ttError = 0;
		double[] ttGradient = new double[paramLength];
		double[][] ttHessian = new double[paramLength][paramLength];

		for (IMoleculeResult result : results) {
			int[] moleculeATs = result.getUpdatedRm().mats;
			int[][] moleculeNPs = result.getUpdatedRm().mnps;
			boolean isDepad = true;

			ttError += result.getTotalError();


			double[] datum = result.getUpdatedRm().datum;

			lla.addData(new ReferenceData(datum[0], result.getHf(),
					ParamGradient.combine(result.getHfDerivs(), info.atomTypes, info.neededParams,
							moleculeATs, moleculeNPs, isDepad),
					ReferenceData.HF_WEIGHT));

			if (datum[1] != 0) {
				lla.addData(new ReferenceData(datum[1], result.getDipole(),
						ParamGradient.combine(result.getDipoleDerivs(), info.atomTypes, info.neededParams,
								moleculeATs, moleculeNPs, isDepad),
						ReferenceData.DIPOLE_WEIGHT));
			}

			if (datum[2] != 0) {
				lla.addData(new ReferenceData(datum[2], result.getIE(),
						ParamGradient.combine(result.getIEDerivs(), info.atomTypes, info.neededParams,
								moleculeATs, moleculeNPs, isDepad),
						ReferenceData.IE_WEIGHT));
			}

			if (result.isExpAvail()) {
				lla.addData(new ReferenceData(0, result.getGeomGradMag(),
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

		logger.info("Total error: {}", ttError);


		// get new search direction
		SimpleMatrix newGradient = new SimpleMatrix(ttGradient);
		SimpleMatrix newHessian = new SimpleMatrix(ttHessian);

//		TROptimizer trOptimizer = new TROptimizer(runInput.lastRunInfo, ttError);
//		lri = trOptimizer.getNewLri();

		IParamOptimizer opt = lla;
		double[] dir = opt.optimize(newHessian, newGradient, null);


		// generating nextInput
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

		InputInfo nextInputInfo = new InputInfo(info.atomTypes, info.neededParams, newNpMap);
		RunnableMolecule[] nextRunRms = new RunnableMolecule[results.length];

		for (int i = 0; i < nextRunRms.length; i++) {
			nextRunRms[i] = results[i].getUpdatedRm();
		}

		RunInput nextInput = new RunInput(nextInputInfo, nextRunRms, lri);


		RunOutput runOutput = new RunOutput(results, sw.getTime(), ttError, ttGradient, ttHessian, runInput,
				nextInput);
		runOutput.finalLambda = opt.getLambda();

		return runOutput;
	}
}
