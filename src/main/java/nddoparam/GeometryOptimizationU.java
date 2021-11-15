package nddoparam;

import org.ejml.simple.SimpleMatrix;

public class GeometryOptimizationU extends GeometryOptimization {

	protected GeometryOptimizationU(SolutionU s) {
		super(s);
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionU) s, i, j);
	}

	protected SimpleMatrix[] findGH() {
		SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine((SolutionU) s);

		SimpleMatrix gradient = matrices[0][0];
		SimpleMatrix hessian;

		// dunno if this is ok for unrestricted
		try {
			 hessian = GeometrySecondDerivative
					.hessianRoutine((SolutionU) s, matrices[1],
							matrices[2]);
		} catch (Exception e) {
			e.printStackTrace();
			hessian = SimpleMatrix.identity(gradient.getNumElements());
		}

		return new SimpleMatrix[]{gradient, hessian};
	}
}
