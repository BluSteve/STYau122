package runcycle.optimize;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.ArrayList;

public class ParamOptimizer {
	private final ArrayList<ReferenceData> datum;
	private double value;
	private static final Logger logger = LogManager.getLogger();

	public ParamOptimizer() {
		this.datum = new ArrayList<>();
	}

	public void addData(ReferenceData data) {
		this.datum.add(data);
		this.value += data.getValue();
	}

	public double[] optimize(SimpleMatrix B, SimpleMatrix gradient) {
		SimpleMatrix searchdir;
		try {
			searchdir = B.pseudoInverse().mult(gradient);
			SimpleMatrix eigenvalues = Utils.symEigen(B)[1];
			logger.info("Hessian eigenvalues: {}", eigenvalues);
		}
		catch (IllegalArgumentException e) {
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
		double lambda = 0;
		double val = 0;
		double[] changes = new double[searchdir.numRows()];

		while (Math.abs(val - value) > 1E-6 && Math.abs(lambda) <= 0.05) {
			lambda += k;
			val = value;
			changes = searchdir.scale(lambda).getDDRM().data;
			value = 0;

			for (ReferenceData d : datum) {
				d.update(changes);
				value += d.getValue();
			}

			if (value > val) {
				k *= -0.5;
			}
		}

		return changes;
	}
}
