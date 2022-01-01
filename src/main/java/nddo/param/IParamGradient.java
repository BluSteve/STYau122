package nddo.param;

import nddo.solution.Solution;

public interface IParamGradient {
	ParamErrorFunction getE();

	Solution getS();

	double[][] getHfDerivs();

	double[][] getDipoleDerivs();

	double[][] getIEDerivs();

	double[][] getGeomDerivs();

	double[][] getTotalGradients();
}
