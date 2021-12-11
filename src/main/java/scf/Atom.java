package scf;

public abstract class Atom {
	protected double[] coordinates;
	protected AtomProperties atomProperties;

	public Atom(AtomProperties atomProperties, double[] coordinates) {
		this.coordinates = coordinates;
		this.atomProperties = atomProperties;
	}

	public double[] getCoordinates() {
		return this.coordinates;
	}

	public AtomProperties getAtomProperties() {
		return this.atomProperties;
	}
}
