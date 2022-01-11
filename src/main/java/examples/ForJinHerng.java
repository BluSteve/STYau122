package examples;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import frontend.FrontendConfig;
import frontend.TxtIO;
import nddo.NDDOAtom;
import nddo.geometry.GeometryOptimization;
import nddo.solution.Solution;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import runcycle.RunIterator;
import runcycle.State;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;
import tools.Batcher;
import tools.Utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class ForJinHerng {
	private static final ScheduledExecutorService progressBar = Executors.newScheduledThreadPool(1);

	public static void main(String[] args) throws IOException, InterruptedException {
		LogManager.getLogger().info("Date compiled: {}", Utils.getResource("version.txt"));
		FrontendConfig.init();

		if (new File("input.txt").isFile() && new File("reference.txt").isFile()) TxtIO.txtToText();

		RunInput input = TxtIO.readInput("molecules.txt");

		Utils.shuffleArray(input.molecules);

		int mLength = input.molecules.length;
		Logger logger = LogManager.getLogger("main");
		logger.info("Molecule count: " + mLength);

		boolean[] isDone = new boolean[mLength];

		AtomicInteger lastLeftCount = new AtomicInteger(mLength);
		AtomicInteger count = new AtomicInteger(0);
		Runnable mLeft = () -> {
			List<String> left = new ArrayList<>(mLength);
			for (int j = 0; j < input.molecules.length; j++) {
				RunnableMolecule molecule = input.molecules[j];
				if (!isDone[molecule.index]) {
					left.add(molecule.debugName());
				}
			}

			int leftCount = left.size();

			if (count.get() % 2 == 1) {
				if (lastLeftCount.get() - leftCount == 0) {
					logger.warn("Stubborn molecule(s) detected, increasing log level...");
					Configurator.setRootLevel(Level.TRACE);
				}
				lastLeftCount.set(leftCount);
			}

			count.incrementAndGet();

			Collections.sort(left);
			logger.info("{}/{} left: {}", leftCount, mLength, left.toString().replace(", ", ","));
		};

		int wait = 30;
		progressBar.scheduleAtFixedRate(mLeft, wait, wait, TimeUnit.SECONDS);

		RunnableMolecule[] updatedRms = new RunnableMolecule[input.molecules.length];
		AtomicInteger doneCount = new AtomicInteger(0);
		Result[] totalResults = Batcher.apply(input.molecules, Result[].class, 1, subset -> {
			Result[] results = new Result[subset.length];

			for (int i = 0; i < subset.length; i++) {
				RunnableMolecule molecule = subset[i];
				try {
					NDDOAtom[] atoms = State.getConverter().convert(molecule.atoms, input.info.npMap);

					Solution s = Solution.of(molecule, atoms);

					GeometryOptimization go = GeometryOptimization.of(s).compute();

					Solution optS = go.getS();

					updatedRms[molecule.index] = new RunnableMolecule(molecule, RunIterator.getAtoms(s),
							molecule.expGeom, molecule.datum);

					double[] numbers = new double[]{s.hf, optS.hf, optS.dipole, -optS.homo};

					results[i] = new Result(molecule.index, molecule.debugName(), numbers);

					isDone[molecule.index] = true;

					molecule.getLogger().info("Finished {}. Initial Hf: {}, final Hf: {}, dipole: {}, IE: {}",
							doneCount.incrementAndGet(), s.hf, optS.hf, optS.dipole, -optS.homo);
				} catch (Exception e) {
					molecule.getLogger().error(molecule.debugName(), e);
					results[i] = new Result(molecule.index, molecule.debugName() + " (errored)", new double[4]);
				}
			}

			return results;
		});

		Arrays.sort(totalResults, Comparator.comparingInt(o -> o.index));

		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		FileWriter writer = new FileWriter("results.json");
		gson.toJson(totalResults, writer);
		writer.close();

		TxtIO.updateMolecules(updatedRms, "updated-molecules.txt");

		if (System.console() != null) {
			System.out.println("Press 'Enter' key to exit.");
			System.console().readLine();
		}
		System.exit(0);
	}

	private static class Result {
		public final int index;
		public final String name;
		public final double[] results;

		public Result(int index, String name, double[] results) {
			this.index = index;
			this.name = name;
			this.results = results;
		}
	}
}
