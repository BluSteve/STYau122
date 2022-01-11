package frontend;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;
import tools.Utils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class Main {
	public static void main(String[] args) throws IOException {
		Logger logger = LogManager.getLogger();
		logger.info("Date compiled: {}", Utils.getResource("version.txt"));
		FrontendConfig.init();

		Files.createDirectories(Path.of("pastinputs"));
		Files.createDirectories(Path.of("outputs"));

		RunInput input;

		try {
			input = JsonIO.readInput("input-override");
			logger.info("Input overriden.");
		} catch (FileNotFoundException e) {
			input = TxtIO.readInput();
		}
		JsonIO.write(input, "original-input");

		RunIterator iterator = new RunIterator(input, FrontendConfig.config.num_runs);

		logger.info("Number of runs = {}", FrontendConfig.config.num_runs);

		int i = 0;
		for (RunOutput ro : iterator) {
			TxtIO.updateInput(ro.nextInput);
			JsonIO.write(ro.nextInput, String.format("pastinputs/%04d-%s", i, ro.nextInput.hash));
			JsonIO.write(ro, String.format("outputs/%04d-%s-%s", i, ro.input.hash, ro.hash));
			i++;
		}

		System.exit(0);
	}
}
