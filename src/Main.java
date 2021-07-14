import nddoparam.NDDOParams;
import nddoparam.mndo.MNDOParams;
import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import runcycle.MoleculeRun;
import runcycle.MoleculeRunR;
import runcycle.MoleculeRunU;
import runcycle.input.InputHandler;
import runcycle.input.RawInput;
import runcycle.input.RawMolecule;
import runcycle.output.MoleculeOutput;
import runcycle.output.OutputHandler;
import scf.AtomHandler;

import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;


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
					"MNDO Parameterization, updated 13 July. " +
							ri.trainingSet +
							" training set (PM7)");

			NDDOParams[] nddoParams =
					new MNDOParams[ri.params.nddoParams.length];
			// TODO change the following line if AM1
			for (int i = 0; i < ri.params.nddoParams.length; i++)
				nddoParams[i] = new MNDOParams(ri.params.nddoParams[i]);

			try {
				List<RawMolecule> requests =
						new ArrayList<>(Arrays.asList(ri.molecules));

				int cores = Runtime.getRuntime().availableProcessors();
				int remainingNonParallel = 5;
				int maxParallel = remainingNonParallel < requests.size() ?
						requests.size() - remainingNonParallel : 1;
				List<RawMolecule> parallelRequests =
						requests.subList(0, maxParallel);
				ForkJoinPool threadPool = new ForkJoinPool(cores);

				List<MoleculeRun> results = threadPool
						.submit(() -> parallelRequests.parallelStream()
								.map(request -> {
									MoleculeRun result = request.restricted ?
											new MoleculeRunR(request,
													nddoParams,
													ri.atomTypes, useHessian) :
											new MoleculeRunU(request,
													nddoParams,
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

				// TODO generalize and put into function for storing
				//  param_length;
				int paramLength = MNDOParams.T1ParamNums.length
						+ (ri.atomTypes.length - 1) *
						MNDOParams.T2ParamNums.length;
				double[] ttGradient =
						new double[paramLength];
				for (MoleculeOutput mo : mos) {
					int k = 0;
					int l = 0;
					for (double[] gradient : mo.gradient.total) {
						int[] paramIndexes;
						if (k == 0) paramIndexes = MNDOParams.T1ParamNums;
						else paramIndexes = MNDOParams.T2ParamNums;
						for (int j : paramIndexes) {
							ttGradient[l] += gradient[j];
							l++;
						}
						k++;
					}
				}
				DoubleMatrix B = findMockHessian(ttGradient,
						ri.params.lastHessian,
						ri.params.lastGradient, ri.params.lastDir,
						paramLength);
				System.out.println(B);

			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		sw.stop();
		System.out.println("Time taken: " + sw.getTime());
	}

	private static DoubleMatrix findMockHessian(double[] newGradient,
												double[] oldHessian,
												double[] oldGradient,
												double[] oldDir, int size) {
		DoubleMatrix s = new DoubleMatrix(oldDir);
		DoubleMatrix hessian = new DoubleMatrix(size, size);
		int count = 1;
		int index = 0;
		while (index < ((size + 1) * size) / 2) {
			for (int i = count - 1; i < size; i++) {
				hessian.put(count - 1, i, oldHessian[index]);
				hessian.put(i, count - 1, oldHessian[index]);
				index++;
			}
			count++;
		}
		DoubleMatrix y =
				new DoubleMatrix(newGradient)
						.sub(new DoubleMatrix(oldGradient));

		double b = y.transpose().mmul(s).get(0);
		DoubleMatrix A = y.mmul(y.transpose()).mmul(1 / b);
		double a = s.transpose().mmul(hessian).mmul(s).get(0);
		DoubleMatrix C = hessian.mmul(s).mmul(s.transpose()).mmul
				(hessian.transpose()).mmul(1 / a);

		return hessian.add(A).sub(C);
	}
}