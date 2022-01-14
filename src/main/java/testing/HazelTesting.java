package testing;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.IExecutorService;
import frontend.FrontendConfig;
import frontend.JsonIO;
import frontend.TxtIO;
import org.apache.logging.log4j.LogManager;
import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;
import runcycle.structs.Serializer;

import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

public class HazelTesting {
	public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
		FrontendConfig.init();

		List<RemoteExecutor> executors = new ArrayList<>();
		String[] ips = {"35.225.69.182"};

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

		for (RemoteExecutor executor : executors) {
			Future<byte[]> future = executor.executorService.submit(new BuildMoleculesTask(pnFile, pFile, mFile));
			RunInput runInput = inflate(future.get(), RunInput.class);
			JsonIO.write(runInput, "remote-input");
			LogManager.getLogger().info("{} {}: {}", executor.coreCount, executor.ip, runInput.molecules.length);


			Future<byte[][]> future2 = executor.executorService.submit(new SolutionTask(runInput));
			byte[][] bytes = future2.get();
			RunOutput ro = inflate(bytes[0], RunOutput.class);
			RunInput nextRunInput = inflate(bytes[1], RunInput.class);
			JsonIO.write(ro, "remote-output");
			JsonIO.write(nextRunInput, "next-run-input");
		}
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

	public static class SolutionTask implements Callable<byte[][]>, Serializable {
		private final byte[] compressedRi;

		public SolutionTask(RunInput ri) {
			this.compressedRi = deflate(ri);
		}

		@Override
		public byte[][] call() {
			RunInput runInput = Serializer.gson.fromJson(Compressor.inflate(compressedRi), RunInput.class);

			RunIterator runIterator = new RunIterator(runInput, FrontendConfig.config.num_runs);

			RunOutput ro = runIterator.next();

			return new byte[][] { deflate(ro), deflate(ro.nextInput)};
		}
	}

	private static byte[] deflate(Object obj) {
		return Compressor.deflate(Serializer.gson.toJson(obj));
	}

	private static <T> T inflate(byte[] bytearr, Class<T> clazz) {
		return Serializer.gson.fromJson(Compressor.inflate(bytearr), clazz);
	}

	public static class BuildMoleculesTask implements Callable<byte[]>, Serializable {
		private final String pnFile, pFile, mFile;

		public BuildMoleculesTask(String pnFile, String pFile, String mFile) {
			this.pnFile = pnFile;
			this.pFile = pFile;
			this.mFile = mFile;
		}

		@Override
		public byte[] call()  {
			return deflate(TxtIO.readInput(pnFile, pFile, mFile));
		}
	}
}
