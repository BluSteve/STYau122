package nddo.param;

import nddo.geometry.GeometryDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionU;

class ParamErrorFunctionU extends ParamErrorFunction {
	protected ParamErrorFunctionU(Solution soln, double refHeat) {
		super(soln, refHeat);
	}

	@Override
	protected double getGradient(int i, int j) {
		return GeometryDerivative.grad((SolutionU) expSoln, i, j);
	}

	@Override
	protected double getDeriv(double[] coeff, int atom2) {
		double deriv = coeff[0] *
				GeometryDerivative.grad((SolutionU) expSoln, atom2, 0)
				+ coeff[1] *
				GeometryDerivative.grad((SolutionU) expSoln, atom2, 1)
				+ coeff[2] *
				GeometryDerivative.grad((SolutionU) expSoln, atom2, 2);
		return 1E-13 * Math.round(deriv * 1E13);
	}
}
