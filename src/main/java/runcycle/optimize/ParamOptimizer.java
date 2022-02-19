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
		newLri.error = error;

		if (lri == null) {
			newLri.trustRadius = 0.01;
		}
		else {
			logger.info("Last run info: {}", lri);
			double changeRatio = (error - lri.error) / lri.expectedChange;
			logger.info("QA validity ratio: {}", changeRatio);
			if (lri.stepSize * 5.0 / 4 > lri.trustRadius && changeRatio >= 0.8) {
				newLri.trustRadius = 1.2 * lri.trustRadius;
			}
			else if (changeRatio <= 0.25) {
				newLri.trustRadius = lri.stepSize / 2;
			}
			else {
				newLri.trustRadius = lri.trustRadius;
			}
		}
	}

	private static double RFOLambda(SimpleMatrix B, SimpleMatrix g) {
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

	private static double lambdamat (SimpleMatrix B, SimpleMatrix g, double alpha) {
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

		SimpleMatrix augmat = SimpleMatrix.identity(B.numRows() + 1);
		augmat.scalei(1 / Math.sqrt(alpha));

		augmat.set(B.numRows(), B.numRows(), 1);

		RFOMat = augmat.mult(RFOMat).mult(augmat);

		return Utils.symEigen(RFOMat)[1].get(0);
	}

	private static double[] RSRFOAlphaLambda(SimpleMatrix B, SimpleMatrix g, double trustradius) {

		double alpha = 1;

		double oldalpha = 0;
		double lambda = 0;
		SimpleMatrix d;
		double ddot = 1;

		while (Math.abs(alpha - oldalpha) > 1E-7 && Math.abs(ddot - trustradius * trustradius) > 1E-10) {
			oldalpha = alpha;
			SimpleMatrix[] matrices = Utils.symEigen(B);
			SimpleMatrix eigenvectors = matrices[0];
			lambda = lambdamat (B, g, alpha);
			d = B.plus(-lambda * alpha, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(g).negativei();
			ddot = d.dot(d);
			double derivative = 0;
			SimpleMatrix dots = g.transpose().mult(eigenvectors);

			for (int i = 0; i < eigenvectors.numCols(); i++) {
				derivative += dots.get(0, i) * dots.get(0, i) * Math.pow(matrices[1].get(i, i) - alpha * lambda, -3);
			}
			derivative *= 2 * lambda / (1 + alpha * d.dot(d));
			logger.info("Microiterations in Progress; RFO shift parameter(s) and validation check(s): {}, {}, {}", lambda, alpha, d.dot(d));

			alpha += 2 * (trustradius * Math.sqrt(d.dot(d)) - d.dot(d)) / derivative;
		}

		if (Math.abs(ddot) < 1E-10) {
			logger.warn("Microiterations did not converge; using standard RFO step (poor)");
			return new double[] {1, RFOLambda(B, g)};
		}


		return new double[] {oldalpha, lambda};
	}

	public double[] optimize(SimpleMatrix B, SimpleMatrix gradient) {
		SimpleMatrix searchdir;
		double[] ls = new double[] {1, 0};
		try {
			//double l = RFOLambda(B, gradient);
			ls = RSRFOAlphaLambda(B, gradient, newLri.trustRadius);

			double l = ls[0] * ls[1];
			logger.info("Aggregated RFO shift parameter: {}", l);
			searchdir = B.plus(-l, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(gradient).negativei();
			SimpleMatrix eigenvalues = Utils.symEigen(B)[1];
			SimpleMatrix eig = eigenvalues.diag().transposei();

			int negCount = 0;
			for (double v : eig.getDDRM().data) if (v < 0) negCount++;

			logger.info("Hessian eigenvalues ({} negative): {}", negCount, eig.toString().strip());
		} catch (IllegalArgumentException e) {
			e.printStackTrace();
			searchdir = gradient.negative();
		}

		double sum = 0;

		for (int i = 0; i < searchdir.numRows(); i++) {
			sum += searchdir.get(i) * searchdir.get(i);
		}

		sum = Math.sqrt(sum);
		searchdir.scalei(1 / sum);

		logger.info("Unscaled step size: {}", sum);


		newLri.stepSize = Math.min(sum, newLri.trustRadius);

		searchdir.scalei(newLri.stepSize);

		double[] changes = searchdir.getDDRM().data;

		logger.info("Final step size: {}", newLri.stepSize);

		newLri.expectedChange = (searchdir.dot(gradient) + 0.5 * searchdir.transpose().mult(B).mult(searchdir).get(0)) / (1 + ls[0] * searchdir.dot(searchdir));

		return changes;
	}

	public LastRunInfo getNewLri() {
		return newLri;
	}
}
