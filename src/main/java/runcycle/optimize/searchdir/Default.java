package runcycle.optimize.searchdir;

import org.ejml.simple.SimpleMatrix;
import runcycle.structs.LastRunInfo;

public class Default implements ISDFinder{
	@Override
	public SimpleMatrix findSD(SimpleMatrix B, SimpleMatrix g, LastRunInfo lri) {
		return B.pseudoInverse().mult(g).negativei();
	}
}
