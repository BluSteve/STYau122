package frontend;

import org.apache.logging.log4j.LogManager;
import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class Main {
	public static void main(String[] args) throws IOException {
		Files.createDirectories(Path.of("pastinputs"));
		Files.createDirectories(Path.of("outputs"));

		RunInput input = TxtIO.readInput();
		JsonIO.write(input, "input");

		RunIterator iterator = new RunIterator(input, FrontendConfig.config.num_runs);

		LogManager.getLogger().info("Number of runs = {}", FrontendConfig.config.num_runs);

		int i = 0;
		for (RunOutput ro : iterator) {
			TxtIO.updateInput(ro.nextInput);
			JsonIO.write(ro.input, String.format("pastinputs/%04d-%s", i, ro.input.hash));
			JsonIO.write(ro, String.format("outputs/%04d-%s-%s", i, ro.input.hash, ro.hash));
			i++;
		}

		System.exit(0);
	}
}
