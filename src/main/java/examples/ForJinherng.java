package examples;

import com.google.gson.Gson;
import frontend.TxtIO;
import nddo.NDDOAtom;
import nddo.geometry.GeometryOptimization;
import nddo.solution.Solution;
import org.apache.logging.log4j.LogManager;
import runcycle.State;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;
import tools.Batcher;
import tools.Utils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.concurrent.atomic.AtomicInteger;

public class ForJinherng {
	private static class Result {
		public final String name;
		public final double[] results;

		public Result(String name, double[] results) {
			this.name = name;
			this.results = results;
		}
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		TxtIO.txtToText();

		RunInput input = TxtIO.readInput();

		Utils.shuffleArray(input.molecules);

		LogManager.getLogger("main").info("Molecule count: " + input.molecules.length);

		AtomicInteger atomicInteger = new AtomicInteger(0);

		Result[] totalResults = new Result[input.molecules.length];

		Batcher.apply(input.molecules, totalResults, subset -> {
			Result[] results = new Result[subset.length];

			for (int i = 0; i < subset.length; i++) {
				RunnableMolecule molecule = subset[i];
				NDDOAtom[] atoms = State.getConverter().convert(molecule.atoms, input.info.npMap);

				Solution s = Solution.of(molecule, atoms);

				GeometryOptimization go = GeometryOptimization.of(s).compute();

				Solution optS = go.getS();

				double[] numbers = new double[]{s.hf, optS.hf, optS.dipole, -optS.homo};

				results[i] = new Result(molecule.debugName(), numbers);

				molecule.getLogger().info("Finished {}. Initial Hf: {}, final Hf: {}, dipole: {}, IE: {}",
						atomicInteger.incrementAndGet(), s.hf, optS.hf, optS.dipole, -optS.homo);
			}

			return results;
		});

		Arrays.sort(totalResults, Comparator.comparing(o -> o.name));

		Gson gson = new Gson();
		FileWriter writer = new FileWriter("results.json");
		gson.toJson(totalResults, writer);
		writer.close();

		if (System.console() != null) {
			System.out.println("Press 'Enter' key to exit.");
			System.console().readLine();
		}
	}
}
