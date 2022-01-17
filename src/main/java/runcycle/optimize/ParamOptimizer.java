package runcycle.optimize;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.ArrayList;

public class ParamOptimizer {
	private static final Logger logger = LogManager.getLogger();
	private final ArrayList<ReferenceData> datum;
	private double value;
	public double lambda;

	private static double lambda(SimpleMatrix h, SimpleMatrix g, int count) {
		if (count == h.numRows()) {
			return 0;
		}

		double initialGuess = h.get(count) - 3;
		double newGuess = initialGuess + 2;
		while (Math.abs(initialGuess - newGuess) > 1E-7) {
			initialGuess = newGuess;
			double f = -initialGuess;
			double fprime = -1;

			for (int i = 0; i < h.numRows(); i++) {
				double v = g.get(i);
				double v1 = h.get(i);
				f += v * v / (initialGuess - v1);
				fprime -= v * v / ((initialGuess - v1) * (initialGuess - v1));
			}

			newGuess = initialGuess - f / fprime;
		}

		if (newGuess != newGuess) {
			throw new IllegalStateException("RFO lambda == null! \n" + h);
		}

		return newGuess;
	}

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
			SimpleMatrix[] ms = Utils.symEigen(B);

			SimpleMatrix h = ms[1].diag();
			SimpleMatrix U = ms[0];

			int counter = 0;
			for (int i = 0; i < h.numRows(); i++) {
				if (Math.abs(h.get(i)) > 1E-5) {
					counter++;
				}
			}

			double l = lambda(h, U.transpose().mult(gradient), h.numRows() - counter);
			searchdir = B.plus(-l, SimpleMatrix.identity(B.numRows())).pseudoInversei().mult(gradient);
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
		double val = 0;
		double[] changes = new double[searchdir.numRows()];

//		while (Math.abs(val - value) > 1E-8) {
//			lambda += k;
//			val = value;
//			changes = searchdir.scale(lambda).getDDRM().data;
//			value = 0;
//
//			for (ReferenceData d : datum) {
//				d.update(changes);
//				value += d.getValue();
//			}
//
//			if (value > val) {
//				k *= -0.5;
//			}
//		}

		logger.info("Final lambda: {}", lambda);

		return changes;
	}
}
