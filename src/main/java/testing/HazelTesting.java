package testing;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.IExecutorService;
import frontend.FrontendConfig;
import frontend.JsonIO;
import frontend.TxtIO;
import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;
import runcycle.structs.Serializer;

import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

public class HazelTesting {
	public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
		FrontendConfig.init();

		ClientConfig clientconf = new ClientConfig();
		clientconf.setClusterName("beJSaHB3AJQVUBa3G7eSptMopuJCRg");
		clientconf.getNetworkConfig().addAddress("34.136.52.126");

		HazelcastInstance h = HazelcastClient.newHazelcastClient(clientconf);
		IExecutorService executorService = h.getExecutorService("serbice");

		List<String> pcsv = Files.readAllLines(Path.of("params.csv"));
		List<String> mtxt = Files.readAllLines(Path.of("inputs/fullch.txt"));

		Future<String> future = executorService.submit(new SolutionTask(pcsv, mtxt));

		RunOutput ro = Serializer.gson.fromJson(future.get(), RunOutput.class);

		JsonIO.write(ro, "remote-output");
	}

	public static class SolutionTask implements Callable<String>, Serializable {
		private final String[] pcsv, mtxt;

		public SolutionTask(List<String> pcsv, List<String> mtxt) {
			this.pcsv = pcsv.toArray(new String[0]);
			this.mtxt = mtxt.toArray(new String[0]);
		}

		@Override
		public String call() throws IOException {
			RunInput runInput = TxtIO.readInput(List.of(pcsv), List.of(mtxt));

			RunIterator runIterator = new RunIterator(runInput, FrontendConfig.config.num_runs);

			RunOutput ro = runIterator.next();

			return Serializer.gson.toJson(ro);
		}
	}
}
