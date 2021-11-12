package nddoparam;

import org.ejml.data.Matrix;
import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;

public class GeometryOptimizationU extends GeometryOptimization {

	protected GeometryOptimizationU(SolutionU s) {
		super(s);
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionU) s, i, j);
	}

	protected SimpleMatrix[] findGH() {
		SimpleMatrix[][] matrices =
				GeometryDerivative.gradientRoutine(s.atoms, (SolutionU) s);

		SimpleMatrix gradient = matrices[0][0];
		SimpleMatrix hessian;

		// dunno if this is ok for unrestricted
		try {
			 hessian = GeometrySecondDerivative
					.hessianRoutine(s.atoms, (SolutionU) s, matrices[1],
							matrices[2]);
		} catch (Exception e) {
			hessian = SimpleMatrix.identity(gradient.getNumElements());
		}

		return new SimpleMatrix[]{gradient, hessian};
	}
}
