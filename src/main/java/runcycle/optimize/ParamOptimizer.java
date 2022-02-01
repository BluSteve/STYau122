package runcycle.optimize;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.structs.LastRunInfo;
import tools.Utils;

public class ParamOptimizer {
	private static final Logger logger = LogManager.getLogger();
	private final LastRunInfo newLri;

	public ParamOptimizer(LastRunInfo lri, double error) {
		newLri = new LastRunInfo();

		if (lri == null) {
			newLri.trustRadius = 0.1;
		}
		else {
			logger.info("Last run info: {}", lri);
			double changeRatio = (error - lri.error) / lri.expectedChange;
			logger.info("QA validity ratio: {}", changeRatio);
			if (lri.stepSize * 5.0 / 4 > lri.trustRadius && changeRatio >= 0.75) {
				newLri.trustRadius = 2 * lri.trustRadius;
			}
			else if (changeRatio <= 0.25) {
				newLri.trustRadius = lri.stepSize / 4;
			}
			else {
				newLri.trustRadius = lri.trustRadius;
			}
		}
	}

	private static double lambda(SimpleMatrix B, SimpleMatrix g) {
		SimpleMatrix RFOMat = new SimpleMatrix(B.numRows() + 1, B.numRows() + 1);

		for (int i = 0; i < B.numRows(); i++) {
			for (int j = i; j < B.numRows(); j++) {
				RFOMat.set(i, j, B.get(i, j));
				RFOMat.set(j, i, B.get(j, i));
			}
		}

		for (int i = 0; i < g.numRows(); i++) {
			RFOMat.set(i, B.numRows(), g.get(i));
			RFOMat.set(B.numRows(), i, g.get(i));
		}

		return Utils.symEigen(RFOMat)[1].get(0);
	}

	public double[] optimize(SimpleMatrix B, SimpleMatrix gradient) {
		SimpleMatrix searchdir;
		try {
			double l = lambda(B, gradient);
			searchdir = B.plus(-l, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(gradient).negativei();
			SimpleMatrix eigenvalues = Utils.symEigen(B)[1];
			if (logger.isInfoEnabled()) {
				SimpleMatrix eig = eigenvalues.diag().transposei();

				int negCount = 0;
				for (double v : eig.getDDRM().data) if (v < 0) negCount++;

				logger.info("Hessian eigenvalues ({} negative): {}", negCount, eig.toString().strip());
			}
		} catch (IllegalArgumentException e) {
			searchdir = gradient.negative();
		}

		double sum = 0;

		for (int i = 0; i < searchdir.numRows(); i++) {
			sum += searchdir.get(i) * searchdir.get(i);
		}

		sum = Math.sqrt(sum);
		searchdir.scalei(1 / sum);

		newLri.stepSize = Math.min(sum, newLri.trustRadius);

		searchdir.scalei(newLri.stepSize);

		double[] changes = searchdir.getDDRM().data;

		logger.info("Final step size: {}", newLri.stepSize);

		newLri.expectedChange = searchdir.dot(gradient) + 0.5 * searchdir.transpose().mult(B).mult(searchdir).get(0);

		return changes;
	}

	public LastRunInfo getNewLri() {
		return newLri;
	}
}
