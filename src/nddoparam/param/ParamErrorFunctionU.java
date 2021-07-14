package nddoparam.param;

import nddoparam.Derivative;
import nddoparam.Solution;
import nddoparam.SolutionU;

public class ParamErrorFunctionU extends ParamErrorFunction {
	public ParamErrorFunctionU(Solution soln, double refHeat) {
		super(soln, refHeat);
	}

	@Override
	protected double getGradient(int i, int j) {
		return Derivative.grad((SolutionU) expSoln, i, j);
	}

	@Override
	protected double getDeriv(double[] coeff, int atom2) {
		double deriv = coeff[0] *
				Derivative.grad((SolutionU) expSoln, atom2, 0)
				+ coeff[1] *
				Derivative.grad((SolutionU) expSoln, atom2, 1)
				+ coeff[2] *
				Derivative.grad((SolutionU) expSoln, atom2, 2);
		return 1E-13 * Math.round(deriv * 1E13);
	}
}
