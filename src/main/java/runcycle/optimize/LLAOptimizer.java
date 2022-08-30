package runcycle.optimize;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.optimize.searchdir.ISDFinder;
import tools.Utils;

import java.util.ArrayList;

public class LLAOptimizer implements IParamOptimizer {
	private static final Logger logger = LogManager.getLogger();
	private final ArrayList<ReferenceData> datum;
	private double valueSum;
	private double lambda;

	public LLAOptimizer() {
		this.datum = new ArrayList<>();
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

	public void addData(ReferenceData data) {
		this.datum.add(data);
		this.valueSum += data.getValue();
	}

	public double[] optimize(SimpleMatrix B, SimpleMatrix gradient, ISDFinder sdFinder) {
		SimpleMatrix searchdir;
		try {
			double l = lambda(B, gradient);
			logger.info("RFO shift parameter: {}", l);

			searchdir = B.plus(-l, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(gradient).negativei();
			SimpleMatrix[] mats = Utils.symEigen(B);

			SimpleMatrix eigenvalues = mats[1];
			SimpleMatrix eig = eigenvalues.diag().transposei();
			int negCount = 0;
			for (double v : eig.getDDRM().data) if (v < 0) negCount++;
			logger.info("Hessian eigenvalues ({} negative): {}", negCount, eig.toString().strip());
		} catch (IllegalArgumentException e) {
			logger.warn("Hessian pinv failed; using gradient instead.");
			searchdir = gradient;
		}

		double sum = 0;

		for (int i = 0; i < searchdir.numRows(); i++) {
			sum += searchdir.get(i) * searchdir.get(i);
		}

		sum = Math.sqrt(sum);
		searchdir = searchdir.scale(1 / sum);


		double k = -0.001;
		lambda = 0;
		double newValueSum = 0;
		double[] changes = new double[searchdir.numRows()];

		while (Math.abs(newValueSum - valueSum) > 1E-8) {
			lambda += k;
			newValueSum = valueSum;
			changes = searchdir.scale(lambda).getDDRM().data;
			valueSum = 0;

			for (ReferenceData d : datum) {
				d.update(changes);
				valueSum += d.getValue();
			}

			if (valueSum > newValueSum) {
				k *= -0.5;
			}
		}

		logger.info("Final lambda: {}", lambda);

		return changes;
	}

	@Override
	public double getLambda() {
		return lambda;
	}
}
