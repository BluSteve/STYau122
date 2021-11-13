package nddoparam;

import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;
import scf.Utils;

public class GeometryOptimizationR extends GeometryOptimization {

	protected GeometryOptimizationR(SolutionR s) {
		super(s);
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionR) s, i, j);
	}

	protected SimpleMatrix[] findGH() {
		DoubleMatrix[][] doubleMatrices =
				GeometryDerivative.gradientRoutine(s.atoms, (SolutionR) s);
		SimpleMatrix[][] matrices = Utils.convertToEJML2D(doubleMatrices);

		SimpleMatrix gradient = matrices[0][0];
		SimpleMatrix hessian;

		try {
			hessian = new SimpleMatrix(GeometrySecondDerivative
					.hessianRoutine(s.atoms, (SolutionR) s, doubleMatrices[1]).toArray2());
		} catch (Exception e) {
			hessian = SimpleMatrix.identity(gradient.getNumElements());
		}

		return new SimpleMatrix[]{gradient, hessian};
	}

}
