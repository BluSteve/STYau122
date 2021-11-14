package nddoparam;

import org.ejml.simple.SimpleMatrix;
import scf.Utils;

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

		double a = 1 / y.transpose().mult(searchdir).get(0);
		double b = searchdir.transpose().mult(B).mult(searchdir).get(0);
		SimpleMatrix m2 = B.mult(searchdir).mult(searchdir.transpose())
				.mult(B.transpose()).scale(b);
		SimpleMatrix m1 = y.mult(y.transpose()).scale(a);

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

		//replaced h summation with simpleh summation

		int counter = 0;
		for (int i = 0; i < h.numRows(); i++) {
			if (Math.abs(h.get(i)) > 1E-5) {
				counter++;
			}
		}

		int numIt = 0;
		while (mag(gradient) > 0.001) {
			// computes new search direction
			double lambda = 0;
			try {
				lambda = lambda(h, U.transpose().mult(gradient),
						h.numRows() - counter);
			} catch (Exception e) {
				e.printStackTrace();
			}
			SimpleMatrix searchDir =
					B.minus(SimpleMatrix.identity(B.numRows()).scale(lambda))
							.invert().mult(gradient).scale(-1);

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

			/*
			 re-compute Hessian if still has not converged after 20 iterations
			*/
			numIt++;
			if (numIt == 15) {
				numIt = 0;
				matrices = findGH();
				B = matrices[1];
			}
			else {
				// creates new gradient
				SimpleMatrix oldGrad = gradient.copy();

				gradient = new SimpleMatrix(s.atoms.length * 3, 1);
				coordIndex = 0;
				for (int a = 0; a < s.atoms.length; a++) {
					for (int i = 0; i < 3; i++) {
						gradient.set(coordIndex, 0, findDerivative(a, i));
						coordIndex++;
					}
				}

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