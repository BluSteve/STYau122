package nddoparam;

import org.jblas.DoubleMatrix;

public class GeometryOptimizationR extends GeometryOptimization {

	public GeometryOptimizationR(NDDOAtom[] atoms, int charge) {
		super(atoms, charge, 1);
	}

	@Override
	protected void updateNDDOSolution() {
		s = new SolutionR(atoms, charge);
	}

	protected double derivative(int i, int j) {
		return Derivative.grad((SolutionR) s, i, j);
	}

	protected DoubleMatrix[] routine() {
		DoubleMatrix[][] matrices =
				Derivative.gradientroutine(atoms, (SolutionR) s);

		DoubleMatrix gradient = matrices[0][0];
		DoubleMatrix hessian;

		try {
			hessian = SecondDerivative
					.hessianroutine(atoms, (SolutionR) s, matrices[1]);
		} catch (Exception e) {
			hessian = DoubleMatrix.eye(gradient.length);
		}

		return new DoubleMatrix[]{gradient, hessian};
	}
}
