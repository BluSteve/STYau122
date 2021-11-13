package nddoparam;

import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;

import static scf.Utils.convertToEJML2D;

public class GeometryOptimizationU extends GeometryOptimization {

	protected GeometryOptimizationU(SolutionU s) {
		super(s);
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionU) s, i, j);
	}

	protected SimpleMatrix[] findGH() {
		DoubleMatrix[][] doubleMatrices =
				GeometryDerivative.gradientRoutine(s.atoms, (SolutionU) s);
		SimpleMatrix[][] matrices = convertToEJML2D(doubleMatrices);

		SimpleMatrix gradient = matrices[0][0];
		SimpleMatrix hessian;

		// dunno if this is ok for unrestricted
		try {
			 hessian = new SimpleMatrix(GeometrySecondDerivative
					.hessianRoutine(s.atoms, (SolutionU) s, doubleMatrices[1],
							doubleMatrices[2]).toArray2());
		} catch (Exception e) {
			hessian = SimpleMatrix.identity(gradient.getNumElements());
		}

		return new SimpleMatrix[]{gradient, hessian};
	}
}
