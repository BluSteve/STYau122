package frontend.json;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import runcycle.structs.RunnableMolecule;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;


public class Main {
	private static final String INPUT_FILENAME = "input";
	private static final int NUM_RUNS = 1;
	private static final int MAX_MOLECULES = 1000;
	private static final boolean isImportLastRun = true;
	private static final Logger logger = LogManager.getLogger();
	private static final ScheduledExecutorService progressBar =
			Executors.newScheduledThreadPool(1);
	private static RawInput ri;
	private static MoleculeOutput[] ranMolecules;

	public static void main(String[] args) throws IOException {
		System.out.println(logger.getLevel());

		if (isImportLastRun) {
			ranMolecules = OutputHandler.importMoleculeOutputs("dynamic-output");
			if (ranMolecules != null) {
				Scanner s = new Scanner(System.in);
				System.out.print(ranMolecules.length + " molecules from last run found, would you like to " +
						"import them? (Y/n) ");
				if (s.next().equals("n")) {
					logger.info("Last ran molecules import cancelled!");
					ranMolecules = null;
				}
				else {
					logger.info("{} molecules successfully imported from last incomplete run!", ranMolecules.length);
				}
				s.close();
			}
			else logger.warn("Last ran molecules import failed!");
		}

		logger.info("NUM_RUNS = " + NUM_RUNS);

		boolean[] isDones = new boolean[MAX_MOLECULES];
		if (logger.isInfoEnabled()) {
			Runnable mLeft = () -> {
				List<String> done = new ArrayList<>();
				List<String> left = new ArrayList<>();

				for (RunnableMolecule rm : ri.molecules) {
					if (isDones[rm.index]) done.add(rm.debugName());
					else left.add(rm.debugName());
				}

				done.sort(String::compareTo);
				left.sort(String::compareTo);

				logger.info("Molecules done ({}): {}", done.size(), done);
				logger.info("Molecules left ({}): {}", left.size(), left);
			};

			int wait = 10;
			progressBar.scheduleAtFixedRate(mLeft, wait, wait,
					TimeUnit.SECONDS);
		}

		StopWatch sw = new StopWatch();
		sw.start();

		sw.stop();
		logger.info("Total time taken: {}", sw.getTime());

		// stop progress bar and terminate program
		progressBar.shutdownNow();
		System.exit(0);
	}
}