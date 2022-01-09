package runcycle.optimize;

public class ReferenceData {
	private final double[] derivatives;
	private final double weight,reference,actual;
	private double value;
	private double updated;
	public static final double HF_WEIGHT = 1;
	public static final double DIPOLE_WEIGHT = 20;
	public static final double GEOM_WEIGHT = 0.5;
	public static final double IE_WEIGHT = 10;

	public ReferenceData(double reference, double actual, double[] derivatives,
						 double weight) {
		this.derivatives = derivatives;
		this.weight = weight;
		this.reference = reference;
		this.actual = actual;
		this.value = (reference - actual) * (reference - actual) * weight * weight;
		this.updated = actual;
	}

	public void update(double[] changes) {
		double change = 0;

		for (int i = 0; i < changes.length; i++) {
			change += derivatives[i] * changes[i];
		}

		this.updated = this.actual + change;

		this.value = (reference - updated) * (reference - updated) * weight * weight;
	}

	public double getValue() {
		return this.value;
	}
}
