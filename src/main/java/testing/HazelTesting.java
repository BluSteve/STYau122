package testing;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.IExecutorService;
import com.sun.management.OperatingSystemMXBean;
import frontend.FrontendConfig;
import frontend.JsonIO;
import frontend.TxtIO;
import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import nddo.structs.MoleculeInfo;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.IMoleculeResult;
import runcycle.RunIterator;
import runcycle.optimize.ParamOptimizer;
import runcycle.optimize.ReferenceData;
import runcycle.structs.*;
import tools.Utils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static runcycle.State.getConverter;

public class HazelTesting {
	private static byte[] deflate(Object obj) {
		return Compressor.deflate(Serializer.gson.toJson(obj));
	}

	private static <T> T inflate(byte[] bytearr, Class<T> clazz) {
		return Serializer.gson.fromJson(Compressor.inflate(bytearr), clazz);
	}

	public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
		FrontendConfig.init();
		Logger logger = LogManager.getLogger();


		// set up Hazelcast
		List<RemoteExecutor> executors = new ArrayList<>();
		String[] ips = {"34.136.5.8", "34.67.122.134", "34.136.23.70"};

		for (String ip : ips) {
			ClientConfig clientconf = new ClientConfig();
			clientconf.setClusterName("beJSaHB3AJQVUBa3G7eSptMopuJCRg");
			clientconf.getNetworkConfig().addAddress(ip + ":5701");
			HazelcastInstance h = HazelcastClient.newHazelcastClient(clientconf);
			IExecutorService executorService = h.getExecutorService("serbice");
			executors.add(new RemoteExecutor(ip, executorService, executorService.submit(new PowerTask()).get()));
		}


		// build initial RunInput object
//		String pnFile = null;
//		String pFile = Files.readString(Path.of("params.csv"));
//		String mFile = Files.readString(Path.of("molecules.txt"));
//
//		RemoteExecutor mainExecutor = executors.get(0);
//		Future<byte[]> future = mainExecutor.executorService.submit(new BuildMoleculesTask(pnFile, pFile, mFile));
//		RunInput runInput = inflate(future.get(), RunInput.class);
//		JsonIO.write(runInput, "remote-input");
//		logger.info("Finished initializing molecules.");

		RunInput runInput = JsonIO.readInput("chnof-run1");

		// creating endingIndices to group molecules by
		int length = runInput.molecules.length;
		double totalpower = 0;
		for (RemoteExecutor executor : executors) totalpower += executor.power;

		List<Integer> endingIndices = new ArrayList<>();
		endingIndices.add(0);
		for (RemoteExecutor executor : executors) {
			int inte = endingIndices.get(endingIndices.size() - 1);
			int end = (int) Math.round(length * (executor.power / totalpower) + inte);
			endingIndices.add(end);
		}
		logger.info("Length: {}, ending indices: {}, powers: {}.", length, endingIndices,
				executors.stream().mapToDouble(e -> e.power));


		// do runs
		int i = FrontendConfig.config.starting_run_num;
		try {
			RunInput currentRunInput = runInput;
			for (; i < FrontendConfig.config.starting_run_num + FrontendConfig.config.num_runs; i++) {
				logger.info("Run number: {}, input hash: {}", i, currentRunInput.hash);

				RunOutput ro = run(executors, endingIndices, currentRunInput);

				logger.info("Run {} time taken: {}, output hash: {}, next input hash: {}\n\n", i, ro.timeTaken,
						ro.hash, ro.nextInputHash);

				currentRunInput = ro.nextInput;

				JsonIO.write(ro, "remote-output");
				JsonIO.write(ro.nextInput, "remote-nextinput");
			}
		} catch (Exception e) {
			logger.error("{} errored!", i, e);
			throw e;
		}
	}

	public static RunOutput run(List<RemoteExecutor> executors, List<Integer> endingIndices, RunInput runInput)
			throws InterruptedException, ExecutionException {
		Logger logger = LogManager.getLogger(runInput.hash);
		RunnableMolecule[] rms = runInput.molecules;
		InputInfo info = runInput.info;
		int nComputers = executors.size();
		StopWatch sw = StopWatch.createStarted();

		// grouping molecules based on coreCount
		Utils.shuffleArray(rms); // shuffles whole rms
		RunnableMolecule[][] rms2d = new RunnableMolecule[nComputers][];
		for (int i = 1; i < endingIndices.size(); i++) {
			RunnableMolecule[] rmsubset = Arrays.copyOfRange(rms, endingIndices.get(i - 1), endingIndices.get(i));
			Arrays.sort(rmsubset, Comparator.comparingInt(rm -> rm.index));
			rms2d[i - 1] = rmsubset;
		}
		Arrays.sort(rms, Comparator.comparingInt(rm -> rm.index)); // deshuffles rms


		// doing the actual computation
		IMoleculeResult[][] results2d = new IMoleculeResult[nComputers][];
		IntStream.range(0, rms2d.length).parallel().forEach(i -> { // multithreaded uploading
			try {
				RemoteExecutor executor = executors.get(i);
				IExecutorService es = executor.executorService;

				RunMoleculesTask task = new RunMoleculesTask(rms2d[i], runInput.info, runInput.hash);

				byte[] resultsBytes = es.submit(task).get();
				results2d[i] = inflate(resultsBytes, IMoleculeResult[].class);

				byte[] logsBytes = es.submit(new LogsTask()).get();
				if (logger.isInfoEnabled()) logger.info("{}:\n{}", executor, Compressor.inflate(logsBytes));
			}
			catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
			}
		});

		IMoleculeResult[] results = Stream.of(results2d).flatMap(Arrays::stream).sorted(
				Comparator.comparingInt(r -> r.getUpdatedRm().index)).toArray(IMoleculeResult[]::new);


		// processing results
		int paramLength = 0; // combined length of all differentiated params
		for (int[] param : info.neededParams) paramLength += param.length;


		// optimizes params based on this run
		ParamOptimizer opt = new ParamOptimizer();
		double ttError = 0;
		double[] ttGradient = new double[paramLength];
		double[][] ttHessian = new double[paramLength][paramLength];

		for (IMoleculeResult result : results) {
			int[] moleculeATs = result.getUpdatedRm().mats;
			int[][] moleculeNPs = result.getUpdatedRm().mnps;
			boolean isDepad = true;

			ttError += result.getTotalError();

			double[] datum = result.getUpdatedRm().datum;

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

		double[] dir = opt.optimize(newHessian, newGradient);


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

		InputInfo nextRunInfo = new InputInfo(info.atomTypes, info.neededParams, newNpMap);
		RunnableMolecule[] nextRunRms = new RunnableMolecule[results.length];

		for (int i = 0; i < nextRunRms.length; i++) {
			nextRunRms[i] = results[i].getUpdatedRm();
		}

		RunInput nextInput = new RunInput(nextRunInfo, nextRunRms);


		return new RunOutput(results, sw.getTime(), ttError, ttGradient, ttHessian, runInput, nextInput);
	}

	public static class RemoteExecutor {
		public final String ip;
		public final IExecutorService executorService;
		public final double power;

		public RemoteExecutor(String ip, IExecutorService executorService, double power) {
			this.ip = ip;
			this.executorService = executorService;
			this.power = power;
		}

		@Override
		public String toString() {
			return "RemoteExecutor{" +
					"ip='" + ip + '\'' +
					", power=" + power +
					'}';
		}
	}

	public static class PowerTask implements Callable<Double>, Serializable {
		@Override
		public Double call() {
			return HazelServer.power;
		}
	}

	public static class BuildMoleculesTask implements Callable<byte[]>, Serializable {
		private final String pnFile, pFile, mFile;

		public BuildMoleculesTask(String pnFile, String pFile, String mFile) {
			this.pnFile = pnFile;
			this.pFile = pFile;
			this.mFile = mFile;
		}

		@Override
		public byte[] call() {
			return deflate(TxtIO.readInput(pnFile, pFile, mFile));
		}
	}

	public static class RunMoleculesTask implements Callable<byte[]>, Serializable {
		private static final Logger logger = LogManager.getLogger();
		private static final OperatingSystemMXBean bean =
				(OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();

		private final byte[] rmsBytes, infoBytes;
		private final String hash;

		public RunMoleculesTask(RunnableMolecule[] rms, InputInfo info, String hash) {
			this.rmsBytes = deflate(rms); // must originally be in sorted order, not necessarily consecutive
			this.infoBytes = deflate(info);
			this.hash = hash;
		}

		@Override
		public byte[] call() {
			// remove past logging
			try {
				new PrintWriter("logs/temp.log").close();
			} catch (FileNotFoundException ignored) {
			}


			// inflate compressed data
			final RunnableMolecule[] rms = inflate(rmsBytes, RunnableMolecule[].class);
			final InputInfo info = inflate(infoBytes, InputInfo.class);


			// start molecule runs
			logger.info("Input hash: {}, nMolecules: {} - Started", hash, rms.length);

			StopWatch sw = StopWatch.createStarted();

			int maxIndex = 0;
			for (RunnableMolecule rm : rms) {
				if (rm.index > maxIndex) maxIndex = rm.index;
			}

			boolean[] isDones = new boolean[maxIndex + 1];
			ScheduledExecutorService progressBar = Executors.newScheduledThreadPool(1);

			if (logger.isInfoEnabled()) {
				AtomicInteger count = new AtomicInteger(0);
				Runnable mLeft = () -> {
					List<RunnableMolecule> left = new ArrayList<>();

					int doneCount = 0;
					for (RunnableMolecule rm : rms) {
						if (!isDones[rm.index]) left.add(rm);
						else doneCount++;
					}
					int leftCount = left.size();

					count.incrementAndGet();

					left.sort(Comparator.comparingInt(rm -> rm.index));

					int totalCount = doneCount + leftCount;
					double progress = 1.0 * doneCount / totalCount;
					long time = sw.getTime(TimeUnit.SECONDS);
					double systemCpuLoad = bean.getSystemCpuLoad();
					double eta = time / progress;
					double percent = 100.0 * progress;

					logger.info("Time: {} s, CPU load: {}, ETA: {} s, {}/{} left ({}% done): {}",
							time,
							String.format("%.2f", systemCpuLoad),
							String.format("%.2f", eta),
							leftCount,
							totalCount,
							String.format("%.2f", percent),
							left.stream().map(MoleculeInfo::debugName).collect(Collectors.toList())
					);
				};

				int wait = FrontendConfig.config.progress_bar_interval;
				progressBar.scheduleAtFixedRate(mLeft, wait, wait, TimeUnit.SECONDS);
			}


			List<RunnableMolecule> rmsList = Arrays.asList(rms);
			Collections.shuffle(rmsList, new Random(Constants.RANDOM_SEED));

			IMoleculeResult[] mResults = rmsList.parallelStream().map(rm -> {
				NDDOAtom[] atoms = getConverter().convert(rm.atoms, info.npMap);
				NDDOAtom[] expGeom = rm.expGeom == null ? null : getConverter().convert(rm.expGeom, info.npMap);

				RunIterator.MoleculeRun mr = new RunIterator.MoleculeRun(rm, atoms, expGeom, rm.datum, true);

				isDones[rm.index] = true;

				return mr;
			}).sorted(Comparator.comparingInt(r -> r.getUpdatedRm().index)).toArray(IMoleculeResult[]::new);

			progressBar.shutdownNow();


			logger.info("Input hash: {}, nMolecules: {} - Finished in {}", hash, rms.length, sw.getTime());
			return deflate(mResults);
		}
	}

	public static class LogsTask implements Callable<byte[]>, Serializable {
		@Override
		public byte[] call() throws IOException {
			String input = Files.readString(Path.of("logs/temp.log")); // todo super bandaid
			return Compressor.deflate(input.replaceAll(String.valueOf(input.charAt(0)), ""));
		}
	}

	@Deprecated
	public static class RunTask implements Callable<byte[][]>, Serializable {
		private final byte[] compressedRi;

		public RunTask(RunInput ri) {
			this.compressedRi = deflate(ri);
		}

		@Override
		public byte[][] call() {
			RunInput runInput = Serializer.gson.fromJson(Compressor.inflate(compressedRi), RunInput.class);

			RunIterator runIterator = new RunIterator(runInput, FrontendConfig.config.num_runs);

			RunOutput ro = runIterator.next();

			return new byte[][]{deflate(ro), deflate(ro.nextInput)};
		}
	}
}
