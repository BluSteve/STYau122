package nddoparam;

import org.jblas.DoubleMatrix;

public class GeometryOptimizationU extends GeometryOptimization {

	public GeometryOptimizationU(NDDOAtom[] atoms, int charge, int mult) {
		super(atoms, charge, mult);
	}

	@Override
	protected void updateNDDOSolution() {
		s = new SolutionU(atoms, charge, mult);
	}

	protected double derivative(int i, int j) {
		return GeometryDerivative.grad((SolutionU) s, i, j);
	}

	protected DoubleMatrix[] routine() {
		DoubleMatrix[][] matrices =
				GeometryDerivative.gradientroutine(atoms, (SolutionU) s);

		DoubleMatrix gradient = matrices[0][0];

		DoubleMatrix hessian = GeometrySecondDerivative
				.hessianroutine(atoms, (SolutionU) s, matrices[1],
						matrices[2]);

		//DoubleMatrix hessian = DoubleMatrix.eye(gradient.length);

		return new DoubleMatrix[]{gradient, hessian};
	}


}
