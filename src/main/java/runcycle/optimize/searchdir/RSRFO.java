package runcycle.optimize.searchdir;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.structs.LastRunInfo;
import tools.Utils;

public class RSRFO implements ISDFinder {
	private static final Logger logger = LogManager.getLogger();

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

			alpha += 2 * (trustradius * Math.sqrt(d.dot(d)) - d.dot(d)) / derivative;
		}

		if (Math.abs(ddot) < 1E-20) {
			logger.warn("Microiterations did not converge; using standard RFO step (poor)");
			return new double[]{1, RFOLambda(B, g)};
		}


		return new double[]{oldalpha, lambda};
	}

	@Override
	public SimpleMatrix findSD(SimpleMatrix B, SimpleMatrix g, LastRunInfo lri) {
		double[] ls = RSRFOAlphaLambda(B, g, lri.trustRadius);
		double l = ls[0] * ls[1]; // RS-RFO step
		return B.plus(-l, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(g).negativei();

	}
}
