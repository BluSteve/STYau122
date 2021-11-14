package nddoparam;

import org.ejml.simple.SimpleMatrix;

public class GeometryOptimizationR extends GeometryOptimization {

	protected GeometryOptimizationR(SolutionR s) {
		super(s);
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionR) s, i, j);
	}

	protected SimpleMatrix[] findGH() {
		SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s.atoms, (SolutionR) s);

		SimpleMatrix gradient = matrices[0][0];
		SimpleMatrix hessian;

		try {
			hessian = GeometrySecondDerivative
					.hessianRoutine(s.atoms, (SolutionR) s, matrices[1]);
		} catch (Exception e) {
			e.printStackTrace();
			hessian = SimpleMatrix.identity(gradient.getNumElements());
		}

		return new SimpleMatrix[]{gradient, hessian};
	}

}
