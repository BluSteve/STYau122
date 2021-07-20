package nddoparam;

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

	protected DoubleMatrix[] findGH() {
		DoubleMatrix[][] matrices =
				GeometryDerivative.gradientRoutine(s.atoms, (SolutionR) s);

		DoubleMatrix gradient = matrices[0][0];
		DoubleMatrix hessian;

		try {
			hessian = GeometrySecondDerivative
					.hessianRoutine(s.atoms, (SolutionR) s, matrices[1]);
		} catch (Exception e) {
			hessian = DoubleMatrix.eye(gradient.length);
		}

		return new DoubleMatrix[]{gradient, hessian};
	}
}
