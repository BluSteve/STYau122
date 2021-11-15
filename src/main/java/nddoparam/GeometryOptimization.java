package nddoparam;

import org.ejml.data.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;

public abstract class GeometryOptimization {
	protected Solution s;

	protected GeometryOptimization(Solution s) {
		this.s = s;
	}

	public static GeometryOptimization of(Solution s) {
		if (s instanceof SolutionR) {
			return new GeometryOptimizationR((SolutionR) s);
		}
		else if (s instanceof SolutionU) {
			return new GeometryOptimizationU((SolutionU) s);
		}
		else throw new IllegalArgumentException(
					"Solution is neither restricted nor unrestricted! " +
							"Molecule: "
							+ s.getRm().index + " " + s.getRm().name);
	}

	private static void logSolution(Solution s) {
		String str = "";
		if (s.getRm() != null) {
			str += s.getRm().index + " " + s.getRm().name + "\n";
		}
		str += "Current heat of formation: " + s.hf + "kcal/mol\n" +
				/*"Current HOMO energy: " + s.homo + " eV\n" +
				"Current energy: " + s.energy + "\n" +*/
				"-----------------------------------------------\n";
		System.out.println(str);
	}

	private static double lambda(SimpleMatrix h, SimpleMatrix g, int count) {
		if (count == h.numRows()) {
			return 0;
		}

		double initialGuess = h.get(count) - 3;
		double newGuess = initialGuess + 2;
		while (Math.abs(initialGuess - newGuess) > 1E-7) {
			initialGuess = newGuess;
			double f = -initialGuess;
			double fprime = -1;

			for (int i = 0; i < h.numRows(); i++) {
				f += g.get(i) * g.get(i) / (initialGuess - h.get(i));
				fprime -= g.get(i) * g.get(i) /
						((initialGuess - h.get(i)) * (initialGuess - h.get(i)));
			}

			newGuess = initialGuess - f / fprime;
		}

		if (newGuess != newGuess) {
			throw new IllegalStateException(
					"RFO lambda == null! \n" + h);
		}

		return newGuess;
	}

	private static SimpleMatrix findNewB(SimpleMatrix B, SimpleMatrix y,
										 SimpleMatrix searchdir) {

		SimpleMatrix yt = y.transpose();
		SimpleMatrix searchdirt = searchdir.transpose();

		double a = 1 / yt.mult(searchdir).get(0);
		double b = searchdirt.mult(B).mult(searchdir).get(0);
		SimpleMatrix m2 = B.mult(searchdir).mult(searchdirt)
				.mult(B.transpose()).scale(b);
		SimpleMatrix m1 = y.mult(yt).scale(a);

		return B.plus(m1).minus(m2);
	}

	private static double mag(SimpleMatrix gradient) {
		double sum = 0;
		for (int i = 0; i < gradient.numRows(); i++) {
			sum += gradient.get(i) * gradient.get(i);
		}

		return Math.sqrt(sum);
	}

	public GeometryOptimization compute() {
		logSolution(s);

		SimpleMatrix[] matrices = findGH();
		SimpleMatrix gradient = matrices[0];
		SimpleMatrix B = matrices[1];

		SimpleMatrix[] ms = Utils.symEigen(B);
		SimpleMatrix h = ms[1].diag();
		SimpleMatrix U = ms[0];

		int counter = 0;
		for (int i = 0; i < h.numRows(); i++) {
			if (Math.abs(h.get(i)) > 1E-5) {
				counter++;
			}
		}

		int numIt = 0;
		while (mag(gradient) > 0.001) {
			// computes new search direction
			double lambda = lambda(h, U.transpose().mult(gradient),
					h.numRows() - counter);

			SimpleMatrix searchDir =
					B.minus(SimpleMatrix.identity(B.numRows()).scale(lambda))
							.invert()
							.mult(gradient)
							.negative();

			if (mag(searchDir) > 0.3) {
				searchDir = searchDir.scale(0.3 / mag(searchDir));
			}

			// updates coordinates of atoms using search direction
			int coordIndex = 0;
			for (NDDOAtom a : s.atoms) {
				for (int i = 0; i < 3; i++) {
					a.getCoordinates()[i] = Math.round(
							(a.getCoordinates()[i] +
									searchDir.get(coordIndex)) *
									1000000000) / 1000000000.0;
					coordIndex++;
				}
			}

			// creates new solution based on updated atom positions
			updateSolution();
			if (s.getRm() != null)
				System.out.println(
						s.getRm().index + " " + s.getRm().name + " Gradient:" +
								" " + mag(gradient));
			logSolution(s);


			// re-compute Hessian if still has not converged after n
			// iterations
			numIt++;
			if (numIt == 15) {
				numIt = 0;
				matrices = findGH();
				B = matrices[1];
			}
			else {
				// creates new gradient
				SimpleMatrix oldGrad = gradient.copy();

				int[][] params = new int[3 * s.atoms.length][];
				double[] results = new double[3 * s.atoms.length];

				coordIndex = 0;
				for (int a = 0; a < s.atoms.length; a++) {
					for (int i = 0; i < 3; i++) {
						params[coordIndex] = new int[]{a, i};
						coordIndex++;
					}
				}

				int elapsedSize = 0;
				double cores = Runtime.getRuntime().availableProcessors();
				int size = Math.max((int) Math.ceil(params.length / cores),
						6);

				List<RecursiveAction> subtasks = new ArrayList<>();
				while (elapsedSize < params.length) {
					int finalElapsedSize = elapsedSize;
					subtasks.add(new RecursiveAction() {
						@Override
						protected void compute() {
							int[][] subset = Arrays.copyOfRange(params,
									finalElapsedSize, Math.min(params.length,
											finalElapsedSize + size));

							double[] output = new double[subset.length];
							for (int i = 0; i < subset.length; i++) {
								output[i] = findDerivative(subset[i][0],
										subset[i][1]);
							}

							System.arraycopy(output, 0, results,
									finalElapsedSize, output.length);
						}
					});
					elapsedSize += size;
				}
				ForkJoinTask.invokeAll(subtasks);

				gradient = new SimpleMatrix(3 * s.atoms.length, 1, true, results);
				SimpleMatrix y = gradient.minus(oldGrad);

				try {
					B = findNewB(B, y, searchDir);
				} catch (Exception e) {
					e.printStackTrace();
					System.err.println("Hessian approximation error!");
					B = SimpleMatrix.identity(s.atoms.length * 3);
				}
			}
			ms = Utils.symEigen(B);

			h = ms[1].diag();
			U = ms[0];
		}
		updateSolution();
		System.out.println("FINAL:");
		logSolution(s);
		return this;
	}

	public Solution getS() {
		return s;
	}

	protected void updateSolution() {
		s = s.withNewAtoms(s.atoms);
	}

	protected abstract SimpleMatrix[] findGH();

	protected abstract double findDerivative(int i, int j);
}