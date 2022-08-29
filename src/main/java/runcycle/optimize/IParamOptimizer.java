package runcycle.optimize;

import org.ejml.simple.SimpleMatrix;
import runcycle.optimize.searchdir.ISDFinder;

public interface IParamOptimizer {
	double[] optimize(SimpleMatrix B, SimpleMatrix g, ISDFinder sdFinder);

	double getLambda();
}
