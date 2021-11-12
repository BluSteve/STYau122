package nddoparam;

import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;
import runcycle.input.RawMolecule;

public class GeometryOptimizationR extends GeometryOptimization {

	protected GeometryOptimizationR(SolutionR s) {
		super(s);
	}

	@Override
	protected void updateSolution() {
		RawMolecule rm = s.getRm();
		s = new SolutionR(s.atoms, s.charge);
		if (rm != null) s.setRm(rm);
	}

	protected double findDerivative(int i, int j) {
		return GeometryDerivative.grad((SolutionR) s, i, j);
	}

	protected SimpleMatrix[] findGH() {
		DoubleMatrix[][] doubleMatrices =
				GeometryDerivative.gradientRoutine(s.atoms, (SolutionR) s);
		SimpleMatrix[][] matrices = convertToEJML2D(doubleMatrices);

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

	private SimpleMatrix[][] convertToEJML2D(DoubleMatrix[][] doubleMatrices) {
		SimpleMatrix[][] matrices = new SimpleMatrix[doubleMatrices.length][];
		for (int i = 0; i < doubleMatrices.length; i++) {
			DoubleMatrix[] dm = doubleMatrices[i];
			matrices[i] = new SimpleMatrix[dm.length];
			for (int j = 0; j < dm.length; j++) {
				matrices[i][j] = new SimpleMatrix(dm[j].toArray2());
			}
		}
		return matrices;
	}
}
