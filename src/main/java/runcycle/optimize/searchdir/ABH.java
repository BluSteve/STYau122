package runcycle.optimize.searchdir;

import org.ejml.simple.SimpleMatrix;
import runcycle.structs.LastRunInfo;
import tools.Utils;

public class ABH implements ISDFinder {
	@Override
	public SimpleMatrix findSD(SimpleMatrix B, SimpleMatrix g, LastRunInfo lri) {
		SimpleMatrix[] mats = Utils.symEigen(B);

		int negCount = 0;
		for (double v : mats[1].diag().transposei().getDDRM().data) if (v < 0) negCount++;

		SimpleMatrix eigenG = mats[0].transpose().mult(g);

		SimpleMatrix searchDir = new SimpleMatrix(mats[0].numRows(), 1);
		for (int i = 0; i < negCount; i++) {
			searchDir.set(i, 0, -eigenG.get(i, 0) / mats[1].get(i, i));
		}
		for (int i = negCount; i < mats[1].numCols(); i++) {
			searchDir.set(i, 0, -eigenG.get(i, 0) / mats[1].get(i, i));
		}

		return searchDir;
	}
}
