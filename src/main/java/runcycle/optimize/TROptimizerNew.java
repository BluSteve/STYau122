package runcycle.optimize;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.optimize.searchdir.ISDFinder;
import runcycle.structs.LastRunInfo;
import tools.Utils;

public class TROptimizerNew implements IParamOptimizer {
	private static final Logger logger = LogManager.getLogger();
	private final LastRunInfo newLri;

	public TROptimizerNew(LastRunInfo lri, double error) {
		newLri = new LastRunInfo();
		newLri.error = error;

		if (lri == null) {
			newLri.trustRadius = 0.01;
		}
		else {
			logger.info("Last run info: {}", lri);
			double changeRatio = (error - lri.error) / lri.expectedChange;
			logger.info("QA validity ratio: {}", changeRatio);
			if (lri.stepSize * 5.0 / 4 > lri.trustRadius && changeRatio >= 0.5 && lri.expectedChange < 5) {
				newLri.trustRadius = 1.1 * lri.trustRadius;
			}
			else if (changeRatio <= 0.25 || lri.expectedChange > 5) {
				newLri.trustRadius = lri.stepSize / 1.2;
			}
			else {
				newLri.trustRadius = lri.trustRadius;
			}
		}
	}

	public double[] optimize(SimpleMatrix B, SimpleMatrix g, ISDFinder sdFinder) {
		SimpleMatrix searchdir;

		try {
			SimpleMatrix eigenvalues = Utils.symEigen(B)[1];
			SimpleMatrix eig = eigenvalues.diag().transposei();

			int negCount = 0;
			for (double v : eig.getDDRM().data) if (v < 0) negCount++;
			logger.info("Hessian eigenvalues ({} negative): {}", negCount, eig.toString().strip());

			searchdir = sdFinder.findSD(B, g, newLri);
		} catch (IllegalArgumentException e) {
			e.printStackTrace();
			searchdir = g.negative();
		}

		newLri.stepSize = newLri.trustRadius;

		searchdir.scalei(newLri.stepSize / searchdir.normF());
		double[] changes = searchdir.getDDRM().data;

		logger.info("Final step size: {}", newLri.stepSize);

		newLri.expectedChange = searchdir.dot(g) + 0.5 * searchdir.transpose().mult(B).mult(searchdir).get(0);

		return changes;
	}

	@Override
	public double getLambda() {
		return newLri.stepSize;
	}

	public LastRunInfo getNewLri() {
		return newLri;
	}
}
