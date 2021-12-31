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
}
