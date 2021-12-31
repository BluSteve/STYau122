package nddo.param;

import nddo.geometry.GeometryDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionR;

class ParamErrorFunctionR extends ParamErrorFunction {
	protected ParamErrorFunctionR(Solution soln, double refHeat) {
		super(soln, refHeat);
	}

	@Override
	protected double getGradient(int i, int j) {
		return GeometryDerivative.grad((SolutionR) expSoln, i, j);
	}
}
