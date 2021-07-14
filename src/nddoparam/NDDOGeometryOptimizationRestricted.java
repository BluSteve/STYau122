package nddoparam;

import org.jblas.DoubleMatrix;

public class NDDOGeometryOptimizationRestricted extends NDDOGeometryOptimization {

	public NDDOGeometryOptimizationRestricted(NDDOAtom[] atoms, int charge) {
		super(atoms, charge, 1);
	}

	@Override
	protected void updateNDDOSolution() {
		s = new NDDOSolutionRestricted(atoms, charge);
	}

	protected double derivative(int i, int j) {
		return NDDODerivative.grad((NDDOSolutionRestricted) s, i, j);
	}

	protected DoubleMatrix[] routine() {
		DoubleMatrix[][] matrices =
				NDDODerivative.gradientroutine(atoms, (NDDOSolutionRestricted) s);

		DoubleMatrix gradient = matrices[0][0];
		DoubleMatrix hessian;

		try {
			hessian = NDDOSecondDerivative
					.hessianroutine(atoms, (NDDOSolutionRestricted) s, matrices[1]);
		} catch (Exception e) {
			hessian = DoubleMatrix.eye(gradient.length);
		}

		return new DoubleMatrix[]{gradient, hessian};
	}
}
