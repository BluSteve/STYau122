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

		List<IExecutorService> executorServices = new ArrayList<>();
		String[] ips = {"localhost", "192.168.31.153"};

		for (String ip : ips) {
			ClientConfig clientconf = new ClientConfig();
			clientconf.setClusterName("beJSaHB3AJQVUBa3G7eSptMopuJCRg");
			clientconf.getNetworkConfig().addAddress(ip);
			clientconf.setProperty("hazelcast.logging.type", "log4j2");
			HazelcastInstance h = HazelcastClient.newHazelcastClient(clientconf);
			IExecutorService executorService = h.getExecutorService("serbice");
			executorServices.add(executorService);
		}

		List<String> pcsv = Files.readAllLines(Path.of("params.csv"));
		List<String> mtxt = Files.readAllLines(Path.of("molecules.txt"));

		for (int i = 0; i < executorServices.size(); i++) {
			Future<String[]> future = executorServices.get(i).submit(new SolutionTask(pcsv, mtxt));
			String[] s = future.get();
			RunOutput ro = Serializer.gson.fromJson(s[0], RunOutput.class);
			RunInput nextInput = Serializer.gson.fromJson(s[1], RunInput.class);
			LogManager.getLogger().info("{}: {}", ips[i], nextInput.info.getParams());
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
