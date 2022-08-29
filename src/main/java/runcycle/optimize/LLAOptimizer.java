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
			SimpleMatrix[] mats = Utils.symEigen(B);

			SimpleMatrix eigenvalues = mats[1];
			SimpleMatrix eig = eigenvalues.diag().transposei();
			int negCount = 0;
			for (double v : eig.getDDRM().data) if (v < 0) negCount++;
			logger.info("Hessian eigenvalues ({} negative): {}", negCount, eig.toString().strip());

			SimpleMatrix grad = mats[0].transpose().mult(gradient);

			searchdir = new SimpleMatrix(mats[0].numRows(), 1);

			for (int i = 0; i < negCount; i++) {
				searchdir.set(i, 0, -grad.get(i, 0) / Math.abs(mats[1].get(i, i)));
			}
			for (int i = negCount; i < mats[1].numCols(); i++) {
				searchdir.set(i, 0, -grad.get(i, 0) / mats[1].get(i, i));
			}

			searchdir = mats[0].mult(searchdir);

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
