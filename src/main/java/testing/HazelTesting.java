package testing;

import com.hazelcast.client.HazelcastClient;
import com.hazelcast.client.config.ClientConfig;
import com.hazelcast.core.HazelcastInstance;
import com.hazelcast.core.IExecutorService;
import frontend.JsonIO;
import frontend.TxtIO;
import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;
import runcycle.structs.Serializer;

import java.io.IOException;
import java.io.Serializable;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

public class HazelTesting {
	public static void main(String[] args) throws IOException, ExecutionException, InterruptedException {
		ClientConfig clientconf = new ClientConfig();
		clientconf.setClusterName("dev");
		clientconf.getNetworkConfig().addAddress("192.168.31.153");

		HazelcastInstance h = HazelcastClient.newHazelcastClient(clientconf);
		IExecutorService executorService = h.getExecutorService("serbice");

		RunInput input = TxtIO.readInput("molecules.txt");

		Future<String> future = executorService.submit(new SolutionTask(input));

		RunOutput ro = Serializer.gson.fromJson(future.get(), RunOutput.class);

		JsonIO.write(ro, "remote-output");
	}

	public static class SolutionTask implements Callable<String>, Serializable {
		private final String inputjson;

		public SolutionTask(RunInput runInput) {
			this.inputjson = Serializer.gson.toJson(runInput);
		}

		@Override
		public String call() {
			RunInput runInput = Serializer.gson.fromJson(inputjson, RunInput.class);

			RunIterator runIterator = new RunIterator(runInput, 1);

			RunOutput ro = runIterator.next();

			return Serializer.gson.toJson(ro);
		}
	}
}
