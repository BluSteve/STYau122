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
import nddo.structs.MoleculeInfo;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import runcycle.IMoleculeResult;
import runcycle.RunIterator;
import runcycle.structs.*;

import java.io.IOException;
import java.io.Serializable;
import java.lang.management.ManagementFactory;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import static runcycle.State.getConverter;

public class HazelTesting {
	public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
		FrontendConfig.init();

		List<RemoteExecutor> executors = new ArrayList<>();
		String[] ips = {"localhost"};

		for (String ip : ips) {
			ClientConfig clientconf = new ClientConfig();
			clientconf.setClusterName("beJSaHB3AJQVUBa3G7eSptMopuJCRg");
			clientconf.getNetworkConfig().addAddress(ip);
			clientconf.setProperty("hazelcast.logging.type", "log4j2");
			HazelcastInstance h = HazelcastClient.newHazelcastClient(clientconf);
			IExecutorService executorService = h.getExecutorService("serbice");
			executors.add(new RemoteExecutor(ip, executorService, executorService.submit(new CoreTask()).get()));
		}

		String pnFile = null;
		String pFile = Files.readString(Path.of("params.csv"));
		String mFile = Files.readString(Path.of("molecules.txt"));

		Logger logger = LogManager.getLogger();
		for (RemoteExecutor executor : executors) {
			Future<byte[]> future = executor.executorService.submit(new BuildMoleculesTask(pnFile, pFile, mFile));
			RunInput runInput = inflate(future.get(), RunInput.class);
			JsonIO.write(runInput, "remote-input");
			logger.info("{} {}: {}", executor.coreCount, executor.ip, runInput.molecules.length);


			Future<byte[]> future2 =
					executor.executorService.submit(new RunMoleculesTask(runInput.molecules, runInput.info));
			byte[] bytes = future2.get();
			IMoleculeResult[] results = inflate(bytes, IMoleculeResult[].class);

			JsonIO.write(results, "remote-results");
			logger.info("run complete: " + bytes.length);
		}
	}

	private static byte[] deflate(Object obj) {
		return Compressor.deflate(Serializer.gson.toJson(obj));
	}

	private static <T> T inflate(byte[] bytearr, Class<T> clazz) {
		return Serializer.gson.fromJson(Compressor.inflate(bytearr), clazz);
	}

	public static class RemoteExecutor {
		public final String ip;
		public final IExecutorService executorService;
		public final int coreCount;

		public RemoteExecutor(String ip, IExecutorService executorService, int coreCount) {
			this.ip = ip;
			this.executorService = executorService;
			this.coreCount = coreCount;
		}
	}

	public static class CoreTask implements Callable<Integer>, Serializable {
		@Override
		public Integer call() {
			return Runtime.getRuntime().availableProcessors();
		}
	}

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

	public static class RunMoleculesTask implements Callable<byte[]>, Serializable {
		private static final Logger logger = LogManager.getLogger();
		private static final OperatingSystemMXBean bean =
				(OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();

		private final byte[] rmsBytes, infoBytes;

		public RunMoleculesTask(RunnableMolecule[] rms, InputInfo info) {
			this.rmsBytes = deflate(rms); // must originally be in sorted order, not necessarily consecutive
			this.infoBytes = deflate(info);
		}

		@Override
		public byte[] call() {
			final RunnableMolecule[] rms = inflate(rmsBytes, RunnableMolecule[].class);
			final InputInfo info = inflate(infoBytes, InputInfo.class);


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
					long time = sw.getTime();
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
							left.stream().map(MoleculeInfo::debugName)
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

			return deflate(mResults);
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
}
