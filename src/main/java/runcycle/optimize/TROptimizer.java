package runcycle.optimize;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.optimize.searchdir.ISDFinder;
import runcycle.structs.LastRunInfo;
import tools.Pow;
import tools.Utils;

public class TROptimizer implements IParamOptimizer {
	private static final Logger logger = LogManager.getLogger();
	private final LastRunInfo newLri;

	public TROptimizer(LastRunInfo lri, double error) {
		newLri = new LastRunInfo();
		newLri.error = error;

		if (lri == null) {
			newLri.trustRadius = 0.1;
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

	private static double lambdamat(SimpleMatrix B, SimpleMatrix g, double alpha) {
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

		while (Math.abs(alpha - oldalpha) > 1E-9 && Math.abs(ddot - trustradius * trustradius) > 1E-15) {
			oldalpha = alpha;
			SimpleMatrix[] matrices = Utils.symEigen(B);
			SimpleMatrix eigenvectors = matrices[0];
			lambda = lambdamat(B, g, alpha);
			d = B.plus(-lambda * alpha, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(g).negativei();
			ddot = d.dot(d);
			double derivative = 0;
			SimpleMatrix dots = g.transpose().mult(eigenvectors);

			for (int i = 0; i < eigenvectors.numCols(); i++) {
				derivative += dots.get(0, i) * dots.get(0, i) * Math.pow(matrices[1].get(i, i) - alpha * lambda, -3);
			}
			derivative *= 2 * lambda / (1 + alpha * d.dot(d));
			//logger.info("Microiterations in Progress; RFO shift parameter(s) and validation check(s): {}, {}, {}",
			// lambda, alpha, d.dot(d));

			alpha += 2 * (trustradius * Math.sqrt(d.dot(d)) - d.dot(d)) / derivative;
		}

		if (Math.abs(ddot) < 1E-20) {
			logger.warn("Microiterations did not converge; using standard RFO step (poor)");
			return new double[]{1, RFOLambda(B, g)};
		}


		return new double[]{oldalpha, lambda};
	}

	private static double lambdaTRMderiv(SimpleMatrix dots, SimpleMatrix eigs, double lambda) {
		double deriv = 0;
		int num = 0;

		if (dots.get(0) == 0) {
			num = 1;
		}

		for (int i = num; i < eigs.numRows(); i++) {
			deriv -= 2 * dots.get(i) * dots.get(i) / Pow.pow(eigs.get(i, i) + lambda, 3);
		}

		return deriv;
	}

	private static double TRMLambda(SimpleMatrix B, SimpleMatrix g, double trustradius) {

		double oldalpha = 0;
		SimpleMatrix d;
		double ddot = 1;

		SimpleMatrix[] mats = Utils.symEigen(B);

		SimpleMatrix dots = g.transpose().mult(mats[0]);

		SimpleMatrix eigs = mats[1];

		double alpha = -eigs.get(0, 0) + 2;

		while (Math.abs(alpha - oldalpha) > 1E-9 && Math.abs(ddot - trustradius * trustradius) > 1E-15) {
			oldalpha = alpha;
			d = B.plus(alpha, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(g).negativei();
			ddot = d.dot(d);
			double derivative = lambdaTRMderiv(dots, eigs, alpha);
			//logger.info("Microiterations in Progress; TRM shift parameter and validation check(s): {}, {}",  alpha,
			// d.dot(d));

			alpha += 2 * (trustradius * Math.sqrt(d.dot(d)) - d.dot(d)) / derivative;
		}

		if (Math.abs(ddot) < 1E-20) {
			logger.warn("Microiterations did not converge; using standard NR step (poor)");
			return 0;
		}


		return -alpha;
	}

	public double[] optimize(SimpleMatrix B, SimpleMatrix gradient, ISDFinder sdFinder) {
		SimpleMatrix searchdir;
		SimpleMatrix oldB = B.copy();

		try {
			SimpleMatrix eigenvalues = Utils.symEigen(B)[1];
			SimpleMatrix eig = eigenvalues.diag().transposei();

			int negCount = 0;
			for (double v : eig.getDDRM().data) if (v < 0) negCount++;
			logger.info("Hessian eigenvalues ({} negative): {}", negCount, eig.toString().strip());

			SimpleMatrix[] mats = Utils.symEigen(B);

			SimpleMatrix grad = mats[0].transpose().mult(gradient);

			searchdir = new SimpleMatrix(mats[0].numRows(), 1);

			for (int i = 0; i < negCount; i++) {
				searchdir.set(i, 0, -grad.get(i, 0) / Math.abs(mats[1].get(i, i)));
			}
			for (int i = negCount; i < mats[1].numCols(); i++) {
				searchdir.set(i, 0, -grad.get(i, 0) / mats[1].get(i, i));
			}

			searchdir = mats[0].mult(searchdir);
//			double[] ls = RSRFOAlphaLambda(oldB, gradient, newLri.trustRadius);
//			double l = ls[0] * ls[1]; //RS-RFO step
//			double l = TRMLambda(oldB, gradient, newLri.trustRadius); //TRM

//			searchdir = B.plus(-l, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(gradient).negativei();
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

		//logger.info("Unscaled step size: {}", sum);


		newLri.stepSize = newLri.trustRadius;


		searchdir.scalei(newLri.stepSize);

		double[] changes = searchdir.getDDRM().data;

		logger.info("Final step size: {}", newLri.stepSize);

		//logger.info ("Gradient values: {}", gradient.toString());


		newLri.expectedChange =
				searchdir.dot(gradient) + 0.5 * searchdir.transpose().mult(oldB).mult(searchdir).get(0);

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
