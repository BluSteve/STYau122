package runcycle.optimize;

import org.ejml.simple.SimpleMatrix;

public interface IParamOptimizer {
	double[] optimize(SimpleMatrix B, SimpleMatrix g);

	double getLambda();
}
