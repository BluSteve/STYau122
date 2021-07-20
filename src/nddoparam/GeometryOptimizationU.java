package nddoparam;

import org.jblas.DoubleMatrix;

public class GeometryOptimizationU extends GeometryOptimization {

	public GeometryOptimizationU(NDDOAtom[] atoms, int charge, int mult) {
		super(atoms, charge, mult);
	}

	@Override
	protected void updateSolution() {
		s = new SolutionU(atoms, charge, mult);
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionU) s, i, j);
	}

	protected DoubleMatrix[] findGH() {
		DoubleMatrix[][] matrices =
				GeometryDerivative.gradientRoutine(atoms, (SolutionU) s);

		DoubleMatrix gradient = matrices[0][0];
		DoubleMatrix hessian;

		// dunno if this is ok for unrestricted
		try {
			 hessian = GeometrySecondDerivative
					.hessianRoutine(atoms, (SolutionU) s, matrices[1],
							matrices[2]);
		} catch (Exception e) {
			hessian = DoubleMatrix.eye(gradient.length);
		}

		return new DoubleMatrix[]{gradient, hessian};
	}
}
