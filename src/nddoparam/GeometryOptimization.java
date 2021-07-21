package nddoparam;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.jblas.Solve;

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
		str += /*"Current heat of formation: " + s.hf + "kcal/mol\n" +
				"Current HOMO energy: " + s.homo + " eV\n" +*/
				"Current energy: " + s.energy + "\n" +
						"-----------------------------------------------\n";
		System.out.println(str);
	}

	private static double lambda(DoubleMatrix h, DoubleMatrix g, int count) {
		if (count == h.length) {
			return 0;
		}

		double initialGuess = h.get(count) - 3;
		double newGuess = initialGuess + 2;
		while (Math.abs(initialGuess - newGuess) > 1E-7) {
			initialGuess = newGuess;
			double f = -initialGuess;
			double fprime = -1;

			for (int i = 0; i < h.rows; i++) {
				f += g.get(i) * g.get(i) / (initialGuess - h.get(i));
				fprime -= g.get(i) * g.get(i) /
						((initialGuess - h.get(i)) * (initialGuess - h.get(i)));
			}

			newGuess = initialGuess - f / fprime;
		}

		if (newGuess != newGuess) {
			throw new IllegalStateException("RFO lambda == null!");
		}

		return newGuess;
	}

	private static DoubleMatrix findNewB(DoubleMatrix B, DoubleMatrix y,
										 DoubleMatrix searchdir) {
		double a = 1 / y.transpose().mmul(searchdir).get(0);
		double b = searchdir.transpose().mmul(B).mmul(searchdir).get(0);
		DoubleMatrix m2 = B.mmul(searchdir).mmul(searchdir.transpose())
				.mmul(B.transpose()).mmul(b);
		DoubleMatrix m1 = y.mmul(y.transpose()).mmul(a);

		return B.add(m1).sub(m2);
	}

	private static double mag(DoubleMatrix gradient) {
		double sum = 0;
		for (int i = 0; i < gradient.length; i++) {
			sum += gradient.get(i) * gradient.get(i);
		}

		return Math.sqrt(sum);
	}

	public GeometryOptimization compute() {
		logSolution(s);

		DoubleMatrix[] matrices = findGH();
		DoubleMatrix gradient = matrices[0];
		DoubleMatrix B = matrices[1];

		DoubleMatrix[] ms = Eigen.symmetricEigenvectors(B);
		DoubleMatrix h = ms[1].diag();
		DoubleMatrix U = ms[0];

		int counter = 0;
		for (int i = 0; i < h.length; i++) {
			if (Math.abs(h.get(i)) > 1E-5) {
				counter++;
			}
		}

		int numIt = 0;
		while (mag(gradient) > 0.001) {
			// computes new search direction
			double lambda = 0;
			try {
				lambda = lambda(h, U.transpose().mmul(gradient),
						h.rows - counter);
			} catch (Exception e) {
				e.printStackTrace();
				System.err.println(h);
			}
			DoubleMatrix searchDir =
					Solve.pinv(B.sub(DoubleMatrix.eye(B.rows).mmul(lambda)))
							.mmul(gradient)
							.mul(-1);
			if (mag(searchDir) > 0.3) {
				searchDir = searchDir.mmul(0.3 / mag(searchDir));
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
			System.out.println(
					s.getRm().index + " " + s.getRm().name + " Gradient: " +
							mag(gradient));
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
				DoubleMatrix oldGrad = gradient.dup();
				gradient = new DoubleMatrix(s.atoms.length * 3, 1);
				coordIndex = 0;
				for (int a = 0; a < s.atoms.length; a++) {
					for (int i = 0; i < 3; i++) {
						gradient.put(coordIndex, 0, findDerivative(a, i));
						coordIndex++;
					}
				}

				// difference of gradients
				DoubleMatrix y = gradient.sub(oldGrad);

				try {
					B = findNewB(B, y, searchDir);
				} catch (Exception e) {
					System.err.println("Hessian approximation error!");
					B = DoubleMatrix.eye(s.atoms.length * 3);
				}
			}

			ms = Eigen.symmetricEigenvectors(B);
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

	protected abstract void updateSolution();

	protected abstract DoubleMatrix[] findGH();

	protected abstract double findDerivative(int i, int j);
}