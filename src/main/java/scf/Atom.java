package scf;

public abstract class Atom {
	protected final double[] coordinates;
	protected final AtomProperties atomProperties;

	public Atom(AtomProperties atomProperties, double[] coordinates) {
		this.coordinates = coordinates.clone();
		this.atomProperties = atomProperties;
	}

	/**
	 * Get coordinates object with no cloning. Can be modified in place!
	 * @return Original coordinates array object.
	 */
	public double[] getCoordinates() {
		return this.coordinates;
	}

	public AtomProperties getAtomProperties() {
		return this.atomProperties;
	}
}
