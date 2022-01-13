package testing;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.IExecutorService;
import frontend.FrontendConfig;
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
		String[] ips = {"localhost", "192.168.31.153"};

		for (String ip : ips) {
			ClientConfig clientconf = new ClientConfig();
			clientconf.setClusterName("beJSaHB3AJQVUBa3G7eSptMopuJCRg");
			clientconf.getNetworkConfig().addAddress(ip);
			clientconf.setProperty("hazelcast.logging.type", "log4j2");
			HazelcastInstance h = HazelcastClient.newHazelcastClient(clientconf);
			IExecutorService executorService = h.getExecutorService("serbice");
			executors.add(new RemoteExecutor(ip, executorService, executorService.submit(new CoreTask()).get()));
		}

		List<String> pcsv = Files.readAllLines(Path.of("params.csv"));
		List<String> mtxt = Files.readAllLines(Path.of("molecules.txt"));

		for (RemoteExecutor executor : executors) {
			Future<String[]> future = executor.executorService.submit(new SolutionTask(pcsv, mtxt));
			String[] s = future.get();
			RunOutput ro = Serializer.gson.fromJson(s[0], RunOutput.class);
			RunInput nextInput = Serializer.gson.fromJson(s[1], RunInput.class);
			LogManager.getLogger().info("{} {}: {}", executor.coreCount, executor.ip, nextInput.info.getParams());
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

	public static class SolutionTask implements Callable<String[]>, Serializable {
		public static final long serialVersionUID = 1234;
		private final String[] pcsv, mtxt;

		public SolutionTask(List<String> pcsv, List<String> mtxt) {
			this.pcsv = pcsv.toArray(new String[0]);
			this.mtxt = mtxt.toArray(new String[0]);
		}

		@Override
		public String[] call() throws IOException {
			RunInput runInput = TxtIO.readInput(List.of(pcsv), List.of(mtxt));

			RunIterator runIterator = new RunIterator(runInput, FrontendConfig.config.num_runs);

			RunOutput ro = runIterator.next();

			return new String[]{Serializer.gson.toJson(ro), Serializer.gson.toJson(ro.nextInput)};
		}
	}
}
