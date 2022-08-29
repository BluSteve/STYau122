package runcycle.optimize.searchdir;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.structs.LastRunInfo;
import tools.Pow;
import tools.Utils;

public class TRM implements ISDFinder {
	private static final Logger logger = LogManager.getLogger();

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

			alpha += 2 * (trustradius * Math.sqrt(d.dot(d)) - d.dot(d)) / derivative;
		}

		if (Math.abs(ddot) < 1E-20) {
			logger.warn("Microiterations did not converge; using standard NR step (poor)");
			return 0;
		}


		return -alpha;
	}

	@Override
	public SimpleMatrix findSD(SimpleMatrix B, SimpleMatrix g, LastRunInfo lri) {
		double l = TRMLambda(B, g, lri.trustRadius);
		return B.plus(-l, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(g).negativei();
	}
}
