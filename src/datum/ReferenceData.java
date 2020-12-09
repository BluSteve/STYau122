package datum;

public abstract class ReferenceData {

    private double[] derivatives;

    private double weight;

    private double reference;

    private double actual;

    private double value;

    private double updated;

    public ReferenceData(double[] derivatives, double weight, double reference, double actual) {

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
