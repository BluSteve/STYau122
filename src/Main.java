import nddoparam.NDDOParams;
import nddoparam.mndo.MNDOParams;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.MoleculeRun;
import runcycle.MoleculeRunR;
import runcycle.MoleculeRunU;
import runcycle.input.InputHandler;
import runcycle.input.RawInput;
import runcycle.input.RawMolecule;
import runcycle.output.MoleculeOutput;
import runcycle.output.OutputHandler;
import scf.AtomHandler;
import scf.Utils;

import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

import static scf.Utils.bohr;

public class Main {

	private static final String INPUT_FILENAME = "input.json";
	private static final String OUTPUT_FILENAME = "output.json";

	public static void main(String[] args) {
		StopWatch sw = new StopWatch();
		sw.start();
//        System.out.close();
		AtomHandler.populateAtoms();

		for (int numRuns = 0; numRuns < 1; numRuns++) {
			boolean useHessian = numRuns % 2 == 0; // Hessian every other run
			InputHandler.processInput(INPUT_FILENAME);
			RawInput ri = InputHandler.ri;
			System.out.println(
					"MNDO Parameterization, updated 13 July. " + ri.trainingSet +
							" training set (PM7)");

			NDDOParams[] nddoParams = new MNDOParams[ri.nddoParams.length];
			// TODO change the following line if AM1
			for (int i = 0; i < ri.nddoParams.length; i++)
				nddoParams[i] = new MNDOParams(ri.nddoParams[i]);

			try {
				List<RawMolecule> requests =
						new ArrayList<>(Arrays.asList(ri.molecules));

				int cores = Runtime.getRuntime().availableProcessors();
				int remainingNonParallel = 5;
				int maxParallel = remainingNonParallel < requests.size() ?
						requests.size() - remainingNonParallel : 1;
				List<RawMolecule> parallelRequests = requests.subList(0, maxParallel);
				ForkJoinPool threadPool = new ForkJoinPool(cores);

				List<MoleculeRun> results = threadPool
						.submit(() -> parallelRequests.parallelStream().map(request -> {
							MoleculeRun result = request.restricted ?
									new MoleculeRunR(request, nddoParams,
											ri.atomTypes, useHessian) :
									new MoleculeRunU(request, nddoParams,
											ri.atomTypes, useHessian);
							return result;
						})).get().collect(Collectors.toList());

				for (RawMolecule request : requests
						.subList(maxParallel, requests.size())) {
					MoleculeRun result = request.restricted ?
							new MoleculeRunR(request, nddoParams, ri.atomTypes,
									useHessian) :
							new MoleculeRunU(request, nddoParams,
									ri.atomTypes,
									useHessian);
					results.add(result);
				}

				MoleculeOutput[] mos = new MoleculeOutput[results.size()];
				for (int i = 0; i < results.size(); i++) {
					mos[i] = OutputHandler.toMoleculeOutput(results.get(i));
				}
				OutputHandler.output(mos, OUTPUT_FILENAME);
				InputHandler.updateInput(ri, INPUT_FILENAME);

			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		sw.stop();
		System.out.println("Time taken: " + sw.getTime());
	}

	private static double[] getHessianUpdateData(Scanner scan) {
		return Utils.toDoubles(scan.nextLine().split(":")[1].split(","));
	}
}