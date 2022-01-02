package nddo.param;

import nddo.solution.Solution;

public interface IParamHessian {
	Solution getS();

	double[][] getHessian();
}
