package nddo.geometry;

import nddo.NDDOAtom;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.apache.logging.log4j.Logger;
import org.ejml.data.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;
import tools.Batcher;
import tools.Utils;

public abstract class GeometryOptimization {
	private final Logger logger;
	protected Solution s;

	protected GeometryOptimization(Solution s) {
		this.s = s;
		logger = s.getRm().getLogger();
	}

	public static GeometryOptimization of(Solution s) {
		if (s instanceof SolutionR) {
			return new GeometryOptimizationR((SolutionR) s);
		}
		else if (s instanceof SolutionU) {
			return new GeometryOptimizationU((SolutionU) s);
		}
		else throw new IllegalArgumentException(
					"Solution is neither restricted nor unrestricted!");
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
			throw new IllegalStateException("RFO lambda == null! \n" + h);
		}

		return newGuess;
	}

	private static SimpleMatrix findNewB(SimpleMatrix B, SimpleMatrix y, SimpleMatrix searchdir) {
		SimpleMatrix yt = y.transpose();
		SimpleMatrix searchdirt = searchdir.transpose();

		double a = 1 / yt.mult(searchdir).get(0);
		double b = searchdirt.mult(B).mult(searchdir).get(0);
		SimpleMatrix m2 = B.mult(searchdir).mult(searchdirt).mult(B.transpose()).scale(b);
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
		logger.debug("initial hf: {}", s.hf);

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
		while (mag(gradient) > 1E-3) {
			// computes new search direction
			double lambda = lambda(h, U.transpose().mult(gradient), h.numRows() - counter);

			SimpleMatrix searchDir = B.minus(SimpleMatrix.identity(B.numRows()).scale(lambda))
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
					a.getCoordinates()[i] = Math.round((a.getCoordinates()[i] + searchDir.get(coordIndex)) *
									1000000000) / 1000000000.0;
					coordIndex++;
				}
			}

			// creates new solution based on updated atom positions
			updateSolution();
			if (logger.isTraceEnabled()) {
				logger.trace("hf: {}, gradient: {}", s.hf, mag(gradient));
			}

			// re-compute Hessian if still has not converged after n
			// iterations
			numIt++;
			if (numIt == 7) {
				numIt = 0;
				matrices = findGH();
				B = matrices[1];
			}
			else {
				// creates new gradient
				SimpleMatrix oldGrad = gradient.copy();

				int[][] params = new int[3 * s.atoms.length][];
				coordIndex = 0;
				for (int a = 0; a < s.atoms.length; a++) {
					for (int i = 0; i < 3; i++) {
						params[coordIndex] = new int[]{a, i};
						coordIndex++;
					}
				}

				double[] results = Batcher.applyDouble(params, subset -> {
					final double[] output = new double[subset.length];
					for (int i = 0; i < subset.length; i++) {
						output[i] = findDerivative(subset[i][0], subset[i][1]);
					}
					return output;
				});

				gradient = new SimpleMatrix(results);

				SimpleMatrix y = gradient.minus(oldGrad);

				try {
					B = findNewB(B, y, searchDir);
				} catch (SingularMatrixException e) {
					s.getRm().getLogger().error("Hessian approximation error!");
					B = SimpleMatrix.identity(s.atoms.length * 3);
				}
			}

			ms = Utils.symEigen(B);
			h = ms[1].diag();
			U = ms[0];
		}

		updateSolution();
		logger.debug("final hf: {}", s.hf);

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