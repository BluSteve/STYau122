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

	public void addData(ReferenceData data) {
		this.datum.add(data);
		this.valueSum += data.getValue();
	}

	public double[] optimize(SimpleMatrix B, SimpleMatrix gradient, ISDFinder sdFinder) {
		SimpleMatrix searchdir;
		try {
			searchdir = B.pseudoInverse().mult(gradient);
			SimpleMatrix eigenvalues = Utils.symEigen(B)[1];
			if (logger.isInfoEnabled()) {
				SimpleMatrix eig = eigenvalues.diag().transposei();

				int negCount = 0;
				for (double v : eig.getDDRM().data) if (v < 0) negCount++;

				logger.info("Hessian eigenvalues ({} negative): {}", negCount, eig.toString().strip());
			}
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

		while (Math.abs(newValueSum - valueSum) > 1E-8 && Math.abs(lambda) <= 5) {
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
