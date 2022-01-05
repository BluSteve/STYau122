package nddo.param;

import nddo.solution.Solution;

public interface IParamHessian {
	Solution getS();

	ParamErrorFunction getE();

	double[][] getHessian();
}
