package runcycle.optimize;

import org.ejml.simple.SimpleMatrix;
import runcycle.IMoleculeResult;
import runcycle.optimize.searchdir.ISDFinder;
import runcycle.structs.InputInfo;
import runcycle.structs.LastRunInfo;
import tools.Utils;

public class TROptimizerNew extends ParamOptimizer {
	public TROptimizerNew(IMoleculeResult[] results, InputInfo info, LastRunInfo lri) {
		super(results, info, lri);

		newLri = new LastRunInfo();
		newLri.error = ttError;

		if (lri == null) {
			newLri.trustRadius = 0.01;
		}
		else {
			logger.info("Last run info: {}", lri);
			double changeRatio = (newLri.error - lri.error) / lri.expectedChange;
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

	public double[] optimize(ISDFinder sdFinder) {
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
	protected void processResult(IMoleculeResult result) {
	}
}
