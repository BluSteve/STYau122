package optimize;

import org.jblas.DoubleMatrix;
import org.jblas.Solve;

import java.util.ArrayList;

public class ParamOptimizer {

	private final ArrayList<ReferenceData> datum;
	public double[] changes;
	private double value;

	public ParamOptimizer() {
		this.datum = new ArrayList<>();
	}

	public void addData(ReferenceData data) {
		this.datum.add(data);
		this.value += data.getValue();
	}

	public double[] optimize(DoubleMatrix B, DoubleMatrix gradient)
			throws Exception {
		DoubleMatrix searchdir = Solve.pinv(B).mmul(gradient);
		double sum = 0;
		for (int i = 0; i < searchdir.rows; i++) {
			sum += searchdir.get(i) * searchdir.get(i);
		}
		sum = Math.sqrt(sum);
		searchdir = searchdir.mmul(1 / sum);

		double k = -0.001;
		double lambda = 0;
		double val = 0;
		this.changes = new double[searchdir.rows];

		int count = 0;
		while (Math.abs(val - value) > 1E-6 && Math.abs(lambda) <= 0.05) {
			count++;
			lambda += k;
			val = value;
			changes = searchdir.dup().mmul(lambda).toArray();
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
