package runcycle.optimize.searchdir;

import org.ejml.simple.SimpleMatrix;
import runcycle.structs.LastRunInfo;

public interface ISDFinder {
	SimpleMatrix findSD(SimpleMatrix B, SimpleMatrix g, LastRunInfo lri);
}
