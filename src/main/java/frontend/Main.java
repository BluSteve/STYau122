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
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Main {
	public static void main(String[] args) throws IOException {
		Logger logger = LogManager.getLogger();
		logger.info("Date compiled: {}", Utils.getResource("version.txt"));
		FrontendConfig.init();

		Files.createDirectories(Path.of("pastinputs"));
		Files.createDirectories(Path.of("outputs"));

		List<RunInput> inputs = new ArrayList<>();
		List<String> filenames = new ArrayList<>();

		try {
			inputs.add(JsonIO.readInput("input-override"));
			filenames.add("overriden-molecules.txt");
			logger.info("Input overriden.");
		} catch (FileNotFoundException e) {
			filenames = Files.walk(Path.of("inputs"))
					.map(path -> path.getFileName().toString())
					.filter(s -> s.endsWith(".txt"))
					.collect(Collectors.toList());


			filenames.forEach(filename -> {
				try {
					logger.info("Input file found: {}", filename);
					inputs.add(TxtIO.readInput(filename));
				} catch (IOException ex) {
					throw new RuntimeException(ex);
				}
			});
		}

		int inputCount = inputs.size();

		logger.info("{} inputs found: {}", inputCount, filenames);

		RunIterator[] iterators = new RunIterator[inputCount];

		for (int i = 0; i < inputCount; i++) {
			iterators[i] = new RunIterator(inputs.get(i), FrontendConfig.config.num_runs);
		}

		String[] filenamesArr = filenames.toArray(new String[0]);

		IntStream.range(0, inputCount).parallel().forEach(index -> {
			try {
				RunIterator iterator = iterators[index];
				String filename = filenamesArr[index];

				int i = FrontendConfig.config.starting_run_num;
				while (iterator.hasNext()) {
					try {
						RunInput current = iterator.getCurrentRunInput();
						JsonIO.write(current, String.format("pastinputs/%04d-%s", i, current.hash));

						RunOutput ro = iterator.next();

						JsonIO.write(ro, String.format("outputs/%04d-%s-%s", i, ro.input.hash, ro.hash));
						TxtIO.updateInput(ro.nextInput, filename);

						i++;
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
				}
			}
			catch (Exception e) {
				logger.error("{} failed!", filenamesArr[index], e);
			}
		});

		System.exit(0);
	}
}
