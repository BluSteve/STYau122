package testing;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.IExecutorService;
import com.sun.management.OperatingSystemMXBean;
import frontend.FrontendConfig;
import frontend.JsonIO;
import frontend.TxtIO;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import nddo.structs.MoleculeInfo;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.IMoleculeResult;
import runcycle.RunIterator;
import runcycle.optimize.ParamOptimizer;
import runcycle.optimize.ReferenceData;
import runcycle.structs.*;
import tools.Utils;

import java.io.*;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static frontend.FrontendConfig.config;
import static runcycle.State.getConverter;

public class HazelTesting {
	private static double[] timeTaken;

	private static byte[] deflateToJson(Object obj) {
		return Compressor.deflate(Serializer.gson.toJson(obj));
	}

	private static <T> T inflateFromJson(byte[] bytearr, Class<T> clazz) {
		return Serializer.gson.fromJson(Compressor.inflate(bytearr), clazz);
	}

	public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
		FrontendConfig.init();
		Logger logger = LogManager.getLogger();

		logger.info("Date compiled: {}", Utils.getResource("version.txt"));

		Files.createDirectories(Path.of("pastinputs"));
		Files.createDirectories(Path.of("outputs"));


		// set up Hazelcast
		List<RemoteExecutor> executors = new ArrayList<>();
		String[] ips = {"104.198.254.246", "localhost",
				"34.122.81.22", "34.71.209.4", "34.75.130.54",
				"35.199.155.191", "35.204.53.185", "35.232.103.247"};

		for (String ip : ips) {
			ClientConfig clientconf = new ClientConfig();
			clientconf.setClusterName("beJSaHB3AJQVUBa3G7eSptMopuJCRg");
			clientconf.getNetworkConfig().addAddress(ip + ":5701");
			HazelcastInstance h = HazelcastClient.newHazelcastClient(clientconf);
			IExecutorService executorService = h.getExecutorService("serbice");
			RemoteExecutor re = new RemoteExecutor(ip, executorService, executorService.submit(new PowerTask()).get());
			executors.add(re);
		}

		executors.sort(Comparator.comparingDouble(e -> -e.power));

		timeTaken = new double[executors.size()];


		// build initial RunInput object
//		String pnFile = Files.readString(Path.of("param-numbers.csv"));
//		String pFile = Files.readString(Path.of("params.csv"));
//		String mFile = Files.readString(Path.of("molecules.txt"));
//
//		RemoteExecutor mainExecutor = executors.get(0);
//		Future<byte[]> future = mainExecutor.executorService.submit(new BuildMoleculesTask(pnFile, pFile, mFile));
//		RunInput runInput = inflateFromJson(future.get(), RunInput.class);
//		JsonIO.write(runInput, "remote-input");
//		logger.info("Finished initializing molecules.");

		RunInput runInput = JsonIO.readInput("remote-input");

		// creating endingIndices to group molecules by
		int length = runInput.molecules.length;
		List<Integer> endingIndices = getEndingIndices(executors, length);

		logger.info("Length: {}, ending indices: {}, executors: {}", length, endingIndices, executors);


		// do runs
		int i = config.starting_run_num;
		int numIt = 1;
		try {
			RunInput currentRunInput = runInput;
			for (; i < config.starting_run_num + config.num_runs; i++) {
				JsonIO.writeAsync(currentRunInput, String.format("pastinputs/%04d-%s", i, currentRunInput.hash));

				logger.info("Run number: {}, input hash: {}", i, currentRunInput.hash);

				RunOutput ro = run(executors, endingIndices, currentRunInput);

				logger.info("Run {} time taken: {}, output hash: {}, next input hash: {}\n\n", i, ro.timeTaken,
						ro.hash, ro.nextInputHash);

				currentRunInput = ro.nextInput;

				JsonIO.writeAsync(ro, String.format("outputs/%04d-%s-%s", i, ro.inputHash, ro.hash));


				System.out.println(Arrays.toString(timeTaken));
				if (numIt % 10 == 0) {
					double max = 0;
					for (double v : timeTaken) if (v > max) max = v;
					for (int j = 0; j < timeTaken.length; j++) timeTaken[j] /= max;

					double spread = Utils.sd(timeTaken);
					if (spread > config.reconf_power_threshold) {
						logger.info("Spread = {}, recalibrating power of machines... (curr={})", spread, endingIndices);

						for (int j = 0; j < executors.size(); j++) {
							executors.get(j).power = 0.5 * executors.get(j).power +
									0.5 * 1 / timeTaken[j] * (endingIndices.get(j + 1) - endingIndices.get(j));
						}
						endingIndices = getEndingIndices(executors, length);

						logger.info("Finished recalibrating power of machines (new={})", endingIndices);

						IntStream.range(0, executors.size()).parallel().forEach(
								j -> executors.get(j).executorService.submit(
										new UpdatePowerTask(executors.get(j).power)));

						logger.info("Uploaded new powers: {}", executors);
					}

					timeTaken = new double[executors.size()];
				}

				numIt++;
			}
		} catch (Exception e) {
			logger.error("{} errored!", i, e);
			throw e;
		}
	}

	public static List<Integer> getEndingIndices(Iterable<RemoteExecutor> executors, int length) {
		double totalpower = 0;
		for (RemoteExecutor executor : executors) totalpower += executor.power;

		List<Integer> endingIndices = new ArrayList<>();
		endingIndices.add(0);
		for (RemoteExecutor executor : executors) {
			int inte = endingIndices.get(endingIndices.size() - 1);
			int end = (int) Math.round(length * (executor.power / totalpower) + inte);
			endingIndices.add(end);
		}
		endingIndices.set(endingIndices.size() - 1, length); // just in case rounding issue

		return endingIndices;
	}

	public static RunOutput run(List<RemoteExecutor> executors, List<Integer> endingIndices, RunInput runInput) {
		Logger logger = LogManager.getLogger(runInput.hash);
		RunnableMolecule[] rms = runInput.molecules;
		InputInfo info = runInput.info;
		int nComputers = executors.size();
		StopWatch sw = StopWatch.createStarted();

		// grouping molecules based on coreCount
		Utils.shuffleArray(rms); // shuffles whole rms, should be the same across runs
		RunnableMolecule[][] rms2d = new RunnableMolecule[nComputers][];
		for (int i = 1; i < endingIndices.size(); i++) {
			RunnableMolecule[] rmsubset = Arrays.copyOfRange(rms, endingIndices.get(i - 1), endingIndices.get(i));
			Arrays.sort(rmsubset, Comparator.comparingInt(rm -> rm.index)); // sorted
			rms2d[i - 1] = rmsubset;
		}
		Arrays.sort(rms, Comparator.comparingInt(rm -> rm.index)); // deshuffles rms


		// doing the actual computation
		Semaphore[] semaphores = new Semaphore[nComputers];
		Queue<IMoleculeResult>[] resultsQs = new Queue[nComputers];
		for (int i = 0; i < nComputers; i++) {
			semaphores[i] = new Semaphore(1);
			resultsQs[i] = new ConcurrentLinkedQueue<>();
		}


		Runnable download = () -> IntStream.range(0, nComputers).parallel().forEach(i -> {
			if (resultsQs[i].size() < rms2d[i].length) {
				try {
					semaphores[i].acquire();

					RemoteExecutor re = executors.get(i);
					List<IMoleculeResult> downloaded = List.of(inflateFromJson(
							re.executorService.submit(new DownloadTask()).get(),
							IMoleculeResult[].class));

					logger.info("Downloaded {} molecules from {}", downloaded.size(), re.ip);

					resultsQs[i].addAll(downloaded);

					semaphores[i].release();
				} catch (InterruptedException | ExecutionException e) {
					e.printStackTrace();
				}
			}
		});

		ScheduledExecutorService ses = Executors.newScheduledThreadPool(1);
		ses.scheduleAtFixedRate(download, 10, 5, TimeUnit.SECONDS);


		IMoleculeResult[][] results2d = new IMoleculeResult[nComputers][];
		AtomicInteger doneCount = new AtomicInteger(0);
		AtomicInteger doneMachineCount = new AtomicInteger(0);
		IntStream.range(0, rms2d.length).parallel().forEach(i -> { // multithreaded uploading
			try {
				RemoteExecutor executor = executors.get(i);
				IExecutorService es = executor.executorService;

				String cachedRmsHash = es.submit(new HashTask()).get();
				String currentRmsHash = Serializer.getHash(rms2d[i]);

				RunMoleculesTask task;
				if (currentRmsHash.equals(cachedRmsHash)) {
					logger.info("{}: using cached runnable molecules.", executor.ip);
					task = new RunMoleculesTask(runInput.info, runInput.hash, executor.ip);
				}
				else {
					logger.info("{}: uploading molecules...", executor.ip);
					task = new RunMoleculesTask(rms2d[i], runInput.info, runInput.hash, executor.ip);
				}


				Pair<Long, byte[]> pair = es.submit(task).get();
				timeTaken[i] += pair.getLeft();
				resultsQs[i].addAll(List.of(inflateFromJson(pair.getRight(), IMoleculeResult[].class)));


				semaphores[i].acquire();
				results2d[i] = resultsQs[i].toArray(new IMoleculeResult[0]);
				semaphores[i].release();
				Arrays.sort(results2d[i], Comparator.comparingInt(r -> r.getUpdatedRm().index));

				if (results2d[i].length != rms2d[i].length) {
					throw new RuntimeException(executor.ip + ": " + rms2d[i].length + " " + results2d[i].length);
				}

				if (logger.isInfoEnabled()) logger.info("{}:\n{}", executor, Compressor.inflate(
						es.submit(new LogsTask()).get()));

				logger.info("{} molecules finished from {}/{} machines", doneCount.addAndGet(results2d[i].length),
						doneMachineCount.incrementAndGet(), nComputers);
			} catch (InterruptedException | ExecutionException e) {
				e.printStackTrace();
			} catch (Exception e) {
				logger.error("", e);
				throw e;
			}
		});

		IMoleculeResult[] results = Stream.of(results2d).flatMap(Arrays::stream).sorted(
				Comparator.comparingInt(r -> r.getUpdatedRm().index)).toArray(IMoleculeResult[]::new);

		ses.shutdownNow();

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


		RunOutput runOutput = new RunOutput(results, sw.getTime(), ttError, ttGradient, ttHessian, runInput,
				nextInput);
		runOutput.finalLambda = opt.lambda;

		return runOutput;
	}

	public static class RemoteExecutor {
		public final String ip;
		public final IExecutorService executorService;
		public double power;

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

	public static class UpdatePowerTask implements Runnable, Serializable {
		private double newPower;

		public UpdatePowerTask(double newPower) {
			this.newPower = newPower;
		}

		@Override
		public void run() {
			try {
				FileWriter pw = new FileWriter("power.txt");
				pw.write(Double.toString(newPower));
				pw.close();
			} catch (IOException e2) {
				e2.printStackTrace();
			}
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
			return deflateToJson(TxtIO.readInput(pnFile, pFile, mFile));
		}
	}

	public static class HashTask implements Callable<String>, Serializable {
		@Override
		public String call() {
			return RunMoleculesTask.cachedRms == null ? "null" : Serializer.getHash(RunMoleculesTask.cachedRms);
		}
	}

	public static class RunMoleculesTask implements Callable<Pair<Long, byte[]>>, Serializable {
		private static final OperatingSystemMXBean bean =
				(OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();
		public static Queue<IMoleculeResult> toDownload = new ConcurrentLinkedQueue<>();
		public static RunnableMolecule[] cachedRms;

		private final byte[] rmsBytes, infoBytes;
		private final String hash, ip;

		public RunMoleculesTask(RunnableMolecule[] rms, InputInfo info, String hash, String ip) {
			this.rmsBytes = deflateToJson(rms); // must originally be in sorted order, not necessarily consecutive
			this.infoBytes = deflateToJson(info);
			this.hash = hash;
			this.ip = ip;
		}

		public RunMoleculesTask(InputInfo info, String hash, String ip) {
			this.rmsBytes = null;
			this.infoBytes = deflateToJson(info);
			this.hash = hash;
			this.ip = ip;
		}

		public synchronized static IMoleculeResult[] emptyDownloadQueue() {
			List<IMoleculeResult> imrList = new LinkedList<>();

			while (!toDownload.isEmpty()) {
				imrList.add(toDownload.remove());
			}

			return imrList.toArray(new IMoleculeResult[0]);
		}

		@Override
		public Pair<Long, byte[]> call() {
			// remove past logging
			try {
				new PrintWriter("logs/temp.log").close();
			} catch (FileNotFoundException ignored) {
			}

			Logger logger = LogManager.getLogger(ip + "_" + hash);


			// inflate compressed data and find whether to use cached values
			final RunnableMolecule[] rms; // not sorted/shuffled in this whole class
			if (rmsBytes == null) {
				logger.info("Using cached runnable molecules.");
				rms = cachedRms;
			}
			else rms = inflateFromJson(rmsBytes, RunnableMolecule[].class);

			final InputInfo info = inflateFromJson(infoBytes, InputInfo.class);


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

				int wait = config.progress_bar_interval;
				progressBar.scheduleAtFixedRate(mLeft, wait, wait, TimeUnit.SECONDS);
			}


			cachedRms = new RunnableMolecule[rms.length];
			IntStream.range(0, rms.length).parallel().forEach(i -> {
				RunnableMolecule rm = rms[i];

				NDDOAtom[] atoms = getConverter().convert(rm.atoms, info.npMap);
				NDDOAtom[] expGeom = rm.expGeom == null ? null : getConverter().convert(rm.expGeom, info.npMap);

				RunIterator.MoleculeRun mr = new RunIterator.MoleculeRun(rm, atoms, expGeom, rm.datum, true);

				isDones[rm.index] = true;

				cachedRms[i] = mr.getUpdatedRm();

				toDownload.add(mr);
			});


			progressBar.shutdownNow();


			logger.info("Input hash: {}, nMolecules: {} - Finished in {}", hash, rms.length, sw.getTime());

			IMoleculeResult[] remaining = emptyDownloadQueue();
			logger.info("Uploading {} remaining molecules...", remaining.length);

			return Pair.of(sw.getTime(), deflateToJson(remaining));
		}
	}

	public static class DownloadTask implements Callable<byte[]>, Serializable {
		@Override
		public byte[] call() {
			IMoleculeResult[] finished = RunMoleculesTask.emptyDownloadQueue();
			LogManager.getLogger().info("Uploading {} molecules...", finished.length);
			return deflateToJson(finished);
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
			this.compressedRi = deflateToJson(ri);
		}

		@Override
		public byte[][] call() {
			RunInput runInput = Serializer.gson.fromJson(Compressor.inflate(compressedRi), RunInput.class);

			RunIterator runIterator = new RunIterator(runInput, config.num_runs);

			RunOutput ro = runIterator.next();

			return new byte[][]{deflateToJson(ro), deflateToJson(ro.nextInput)};
		}
	}
}
