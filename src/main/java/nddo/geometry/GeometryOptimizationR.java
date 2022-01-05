package nddo.geometry;

import nddo.solution.SolutionR;
import org.ejml.data.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;

public class GeometryOptimizationR extends GeometryOptimization {
	protected GeometryOptimizationR(SolutionR s) {
		super(s);
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionR) s, i, j);
	}

	protected SimpleMatrix[] findGH() {
		SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine((SolutionR) s);

		SimpleMatrix gradient = matrices[0][0];
		SimpleMatrix hessian;

		try {
			hessian = GeometrySecondDerivative.hessianRoutine((SolutionR) s, matrices[1]);
		} catch (SingularMatrixException | IllegalStateException e) {
			s.rm.getLogger().warn("hessianRoutine (R) failed: {}", e.toString());
			hessian = SimpleMatrix.identity(gradient.getNumElements());
		}

		return new SimpleMatrix[]{gradient, hessian};
	}
}
