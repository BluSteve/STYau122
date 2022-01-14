package testing;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.IExecutorService;
import frontend.FrontendConfig;
import frontend.JsonIO;
import frontend.TxtIO;
import org.apache.logging.log4j.LogManager;
import runcycle.structs.RunInput;
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
			Future<String> future = executor.executorService.submit(new BuildMoleculesTask(pnFile, pFile, mFile));
			String s = future.get();
			RunInput runInput = Serializer.gson.fromJson(s, RunInput.class); // not deserializing correctly
			JsonIO.write(runInput, "remote-input");
			LogManager.getLogger().info("{} {}: {}", executor.coreCount, executor.ip, runInput.molecules.length);
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

//	public static class SolutionTask implements Callable<String[]>, Serializable {
//		private final String riJson;
//
//		public SolutionTask(RunInput ri) {
//			this.rmsJson = Serializer.gson.toJson(rms);
//		}
//
//		@Override
//		public String[] call() throws IOException {
//			RunInput runInput = TxtIO.readInput(List.of(pcsv), List.of(mtxt));
//
//			RunIterator runIterator = new RunIterator(runInput, FrontendConfig.config.num_runs);
//
//			RunOutput ro = runIterator.next();
//
//			return new String[]{Serializer.gson.toJson(ro), Serializer.gson.toJson(ro.nextInput)};
//		}
//	}

	public static class BuildMoleculesTask implements Callable<String>, Serializable {
		private final String pnFile, pFile, mFile;

		public BuildMoleculesTask(String pnFile, String pFile, String mFile) {
			this.pnFile = pnFile;
			this.pFile = pFile;
			this.mFile = mFile;
		}

		@Override
		public String call()  {
			return Serializer.gson.toJson(TxtIO.readInput(pnFile, pFile, mFile));
		}
	}
}
