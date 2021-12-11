package scf;

public abstract class Atom {
	protected final double[] coordinates;
	protected final AtomProperties atomProperties;

	/**
	 * Coordinates is like passing in 3 double values.
	 */
	public Atom(AtomProperties atomProperties, double[] coordinates) {
		this.coordinates = coordinates.clone();
		this.atomProperties = atomProperties;
	}

	public double[] getCoordinates() {
		return this.coordinates;
	}

	public AtomProperties getAtomProperties() {
		return this.atomProperties;
	}
}
