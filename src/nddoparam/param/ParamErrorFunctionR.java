package nddoparam.param;

import nddoparam.GeometryDerivative;
import nddoparam.Solution;
import nddoparam.SolutionR;

public class ParamErrorFunctionR extends ParamErrorFunction {
	public ParamErrorFunctionR(Solution soln, double refHeat) {
		super(soln, refHeat);
	}

	@Override
	protected double getGradient(int i, int j) {
		return GeometryDerivative.grad((SolutionR) expSoln, i, j);
	}

	@Override
	protected double getDeriv(double[] coeff, int atom2) {
		double deriv = coeff[0] *
				GeometryDerivative.grad((SolutionR) expSoln, atom2, 0)
				+ coeff[1] *
				GeometryDerivative.grad((SolutionR) expSoln, atom2, 1)
				+ coeff[2] *
				GeometryDerivative.grad((SolutionR) expSoln, atom2, 2);
		return 1E-13 * Math.round(deriv * 1E13);
	}
}
