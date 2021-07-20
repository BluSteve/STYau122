package nddoparam;

import org.jblas.DoubleMatrix;

public class GeometryOptimizationR extends GeometryOptimization {

	public GeometryOptimizationR(NDDOAtom[] atoms, int charge) {
		super(atoms, charge, 1);
	}

	@Override
	protected void updateSolution() {
		s = new SolutionR(atoms, charge);
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionR) s, i, j);
	}

	protected DoubleMatrix[] findGH() {
		DoubleMatrix[][] matrices =
				GeometryDerivative.gradientRoutine(atoms, (SolutionR) s);

		DoubleMatrix gradient = matrices[0][0];
		DoubleMatrix hessian;

		try {
			hessian = GeometrySecondDerivative
					.hessianRoutine(atoms, (SolutionR) s, matrices[1]);
		} catch (Exception e) {
			hessian = DoubleMatrix.eye(gradient.length);
		}

		return new DoubleMatrix[]{gradient, hessian};
	}
}
