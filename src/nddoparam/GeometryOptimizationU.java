package nddoparam;

import org.jblas.DoubleMatrix;

public class GeometryOptimizationU extends GeometryOptimization {

	protected GeometryOptimizationU(SolutionU s) {
		super(s);
	}

	@Override
	protected void updateSolution() {
		s.compute();
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionU) s, i, j);
	}

	protected DoubleMatrix[] findGH() {
		DoubleMatrix[][] matrices =
				GeometryDerivative.gradientRoutine(s.atoms, (SolutionU) s);

		DoubleMatrix gradient = matrices[0][0];
		DoubleMatrix hessian;

		// dunno if this is ok for unrestricted
		try {
			 hessian = GeometrySecondDerivative
					.hessianRoutine(s.atoms, (SolutionU) s, matrices[1],
							matrices[2]);
		} catch (Exception e) {
			hessian = DoubleMatrix.eye(gradient.length);
		}

		return new DoubleMatrix[]{gradient, hessian};
	}
}
