package nddoparam;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.jblas.Solve;

public abstract class GeometryOptimization {
	public Solution s;
	protected int charge;
	protected NDDOAtom[] atoms;
	protected int mult;
	protected DoubleMatrix gradient;

	public GeometryOptimization(NDDOAtom[] atoms, int charge, int mult) {
		this.atoms = atoms;
		this.charge = charge;
		this.mult = mult;
		updateSolution();

		logSolution(s);

		DoubleMatrix[] matrices = findGH();

		gradient = matrices[0];
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
			System.out.println("Gradient: " + mag(gradient));

			// computes new search directions
			double lambda =
					lambda(h, U.transpose().mmul(gradient), h.rows - counter);
			DoubleMatrix searchDir =
					Solve.pinv(B.sub(DoubleMatrix.eye(B.rows).mmul(lambda)))
							.mmul(gradient)
							.mul(-1);
			if (mag(searchDir) > 0.3) {
				searchDir = searchDir.mmul(0.3 / mag(searchDir));
			}

			// updates coordinates of atoms
			int coordIndex = 0;
			for (NDDOAtom a : atoms) {
				for (int i = 0; i < 3; i++) {
					a.getCoordinates()[i] = Math.round(
							(a.getCoordinates()[i] +
									searchDir.get(coordIndex)) *
									1000000000) / 1000000000.0;
					coordIndex++;
				}
			}

			updateSolution();
			logSolution(s);

			// Re-compute Hessian if still has not converged
			numIt++;
			if (numIt == 20) {
				numIt = 0;
				matrices = findGH();
				B = matrices[1];
			}
			else {
				// creates new gradient
				DoubleMatrix oldGrad = gradient.dup();
				gradient = new DoubleMatrix(atoms.length * 3, 1);
				coordIndex = 0;
				for (int a = 0; a < atoms.length; a++) {
					for (int i = 0; i < 3; i++) {
						gradient.put(coordIndex, 0, findDerivative(a, i));
						coordIndex++;
					}
				}

				// difference of gradients
				DoubleMatrix y = gradient.sub(oldGrad);

				try {
					B = findB(B, y, searchDir);
				} catch (Exception e) {
					System.err.println("Hessian approximation error!");
					B = DoubleMatrix.eye(atoms.length * 3);
				}
			}

			ms = Eigen.symmetricEigenvectors(B);
			h = ms[1].diag();
			U = ms[0];
		}
		updateSolution();
		System.out.println("FINAL:");
		logSolution(s);
	}

	private static void logSolution(Solution s) {
		System.out.println(s.getRm().index + " " + s.getRm().name + "\n" +
				"Current heat of formation: " + s.hf + "kcal/mol\n" +
				"Current HOMO energy: " + s.homo + " eV\n" +
				"Current energy: " + s.energy + "\n" +
				"-----------------------------------------------\n");
	}

	private static double lambda(DoubleMatrix h, DoubleMatrix g, int count) {
		if (count == h.length) {
			return 0;
		}

		double initialguess = h.get(count) - 3;
		double newguess = initialguess + 2;
		while (Math.abs(initialguess - newguess) > 1E-7) {
			initialguess = newguess;
			double f = -initialguess;
			double fprime = -1;

			for (int i = 0; i < h.rows; i++) {
				f += g.get(i) * g.get(i) / (initialguess - h.get(i));
				fprime -= g.get(i) * g.get(i) /
						((initialguess - h.get(i)) * (initialguess - h.get(i)));
			}

			newguess = initialguess - f / fprime;
		}

		return newguess;
	}

	private static DoubleMatrix findB(DoubleMatrix B, DoubleMatrix y,
									  DoubleMatrix searchdir) {
		double a = 1 / y.transpose().mmul(searchdir).get(0);
		double b = searchdir.transpose().mmul(B).mmul(searchdir).get(0);
		DoubleMatrix m2 =
				B.mmul(searchdir).mmul(searchdir.transpose())
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

	public NDDOAtom[] getAtoms() {
		return this.atoms;
	}

	protected abstract void updateSolution();

	protected abstract DoubleMatrix[] findGH();

	protected abstract double findDerivative(int i, int j);
}