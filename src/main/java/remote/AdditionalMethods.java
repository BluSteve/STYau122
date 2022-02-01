package remote;

import com.sun.management.OperatingSystemMXBean;
import frontend.FrontendConfig;
import frontend.TxtIO;
import nddo.NDDOAtom;
import nddo.structs.MoleculeInfo;
import node.Node;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import runcycle.IMoleculeResult;
import runcycle.RunIterator;
import runcycle.structs.InputInfo;
import runcycle.structs.RunnableMolecule;
import runcycle.structs.Serializer;
import shared.Constants;
import shared.Message;
import shared.MethodContainer;
import tools.Byter;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static frontend.FrontendConfig.config;
import static remote.Utils.fromJsonBytes;
import static remote.Utils.toJsonBytes;
import static runcycle.State.getConverter;

@SuppressWarnings("unused")
public class AdditionalMethods extends MethodContainer {
	private static final OperatingSystemMXBean bean =
			(OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();
	private static final Logger logger = LogManager.getLogger();
	private static final int DOWNLOAD_THRESHOLD = 20;
	private final Queue<IMoleculeResult> toDownload = new ConcurrentLinkedQueue<>();
	private RunnableMolecule[] cachedRms;

	public AdditionalMethods(Node node) {
		super(node);
	}

	private static void powBench() {
		IntStream.range(0, 1_000_000).parallel().forEach(i -> {
			Random r = new Random();
			double product = 1;
			for (int j = 0; j < 1000; j++) {
				product *= Math.pow(r.nextGaussian(), r.nextGaussian());
			}
		});
	}

	public byte[] ping() {
		return "pong!".getBytes(StandardCharsets.UTF_8);
	}

	public byte[] benchmark() {
		final String powertxt = "power.txt";
		double power;

		try {
			if (Files.exists(Path.of(powertxt))) {
				power = Double.parseDouble(Files.readString(Path.of(powertxt)));
			}
			else {
				powBench();

				long time = System.nanoTime();

				powBench();

				power = 10 / ((System.nanoTime() - time) / 1e9);

				FileWriter pw = new FileWriter(powertxt);
				pw.write(Double.toString(power));
				pw.close();
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}


		return Byter.toBytes(power);
	}

	public void updatePower(byte[] bytes) {
		double newPower = Byter.toDouble(bytes);

		try {
			FileWriter pw = new FileWriter("power.txt");
			pw.write(Double.toString(newPower));
			pw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public byte[] buildMolecules(byte[] bytes) {
		AdvancedMachine.MoleculeInputFiles inputFiles = fromJsonBytes(bytes, AdvancedMachine.MoleculeInputFiles.class);

		return Utils.toJsonBytes(TxtIO.readInput(inputFiles.pnFile, inputFiles.pFile, inputFiles.mFile));
	}

	public byte[] uploadMolecules() {
		IMoleculeResult[] finished = emptyDownloadQueue();
		logger.info("Uploading {} molecules...", finished.length);
		return Utils.toJsonBytes(finished);
	}

	public byte[] getLogs() throws IOException {
		String input = Files.readString(Path.of("logs/temp.log")); // todo super bandaid
		return input.replaceAll(String.valueOf(input.charAt(0)), "").getBytes(StandardCharsets.UTF_8);
	}

	public byte[] getMoleculesHash() {
		String hash = cachedRms == null ? "null" : Serializer.getHash(cachedRms);

		return hash.getBytes(StandardCharsets.UTF_8);
	}

	public byte[] runMolecules(byte[] bytes) {
		FrontendConfig.init();
		AdvancedMachine.MoleculesSubset subset = fromJsonBytes(bytes, AdvancedMachine.MoleculesSubset.class);

		// remove past logging
		try {
			new PrintWriter("logs/temp.log").close();
		} catch (FileNotFoundException ignored) {
		}

		Logger logger = LogManager.getLogger(subset.ip + "_" + subset.inputHash);


		// inflate compressed data and find whether to use cached values
		final RunnableMolecule[] rms; // not sorted/shuffled in this whole class
		if (subset.rms == null) {
			logger.info("Using cached runnable molecules.");
			rms = cachedRms;
		}
		else rms = subset.rms;

		final InputInfo info = subset.info;


		StopWatch sw = StopWatch.createStarted();
		logger.info("Input hash: {}, nMolecules: {} - Started", subset.inputHash, rms.length);


		// progress bar
		int maxIndex = 0;
		for (RunnableMolecule rm : rms) {
			if (rm.index > maxIndex) maxIndex = rm.index;
		}

		boolean[] isDones = new boolean[maxIndex + 1];
		ScheduledExecutorService progressBar = Executors.newScheduledThreadPool(1);

		if (logger.isInfoEnabled()) {
			AtomicInteger count = new AtomicInteger(0);
			Runnable mLeft = () -> {
				List<RunnableMolecule> left = new ArrayList<>();

				int doneCount = 0;
				for (RunnableMolecule rm : rms) {
					if (!isDones[rm.index]) left.add(rm);
					else doneCount++;
				}
				int leftCount = left.size();

				count.incrementAndGet();

				left.sort(Comparator.comparingInt(rm -> rm.index));

				int totalCount = doneCount + leftCount;
				double progress = 1.0 * doneCount / totalCount;
				long time = sw.getTime(TimeUnit.SECONDS);
				double systemCpuLoad = bean.getSystemCpuLoad();
				double eta = time / progress;
				double percent = 100.0 * progress;

				logger.info("Time: {} s, CPU load: {}, ETA: {} s, {}/{} left ({}% done): {}",
						time,
						String.format("%.2f", systemCpuLoad),
						String.format("%.2f", eta),
						leftCount,
						totalCount,
						String.format("%.2f", percent),
						left.stream().map(MoleculeInfo::debugName).collect(Collectors.toList())
				);
			};

			int wait = config.progress_bar_interval;
			progressBar.scheduleAtFixedRate(mLeft, wait, wait, TimeUnit.SECONDS);
		}


		// dynamic download
		AtomicReference<Exception> error = new AtomicReference<>();
		ExecutorService downloadService = Executors.newFixedThreadPool(1);
		downloadService.submit(() -> {
			while (!Thread.currentThread().isInterrupted()) {
				if (toDownload.size() >= DOWNLOAD_THRESHOLD) {
					try {
						new Message(Constants.response, uploadMolecules()).writeTo(node.out);
					} catch (IOException e) {
						error.set(e);
						break;
					}
				}
			}
		});


		// start molecule runs
		cachedRms = new RunnableMolecule[rms.length];
		IntStream.range(0, rms.length).parallel().forEach(i -> {
			Exception e = error.get();
			if (e == null) {
				RunnableMolecule rm = rms[i];

				NDDOAtom[] atoms = getConverter().convert(rm.atoms, info.npMap);
				NDDOAtom[] expGeom = rm.expGeom == null ? null : getConverter().convert(rm.expGeom, info.npMap);

				RunIterator.MoleculeRun mr = new RunIterator.MoleculeRun(rm, atoms, expGeom, rm.datum, true);

				isDones[rm.index] = true;

				cachedRms[i] = mr.getUpdatedRm();

				toDownload.add(mr);
			}
			else throw new RuntimeException(e);
		});


		progressBar.shutdownNow();
		downloadService.shutdownNow();


		logger.info("Input hash: {}, nMolecules: {} - Finished in {}", subset.inputHash, rms.length, sw.getTime());

		IMoleculeResult[] remaining = emptyDownloadQueue();
		logger.info("Uploading {} remaining molecules...", remaining.length);

		return toJsonBytes(new AdvancedMachine.SubsetResult(sw.getTime(), remaining));
	}

	private synchronized IMoleculeResult[] emptyDownloadQueue() {
		List<IMoleculeResult> imrList = new LinkedList<>();

		while (!toDownload.isEmpty()) {
			imrList.add(toDownload.remove());
		}

		return imrList.toArray(new IMoleculeResult[0]);
	}
}
