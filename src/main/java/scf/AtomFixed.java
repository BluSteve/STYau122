package scf;
// AtomProperties <--HAS-A-- AtomFixed <--IS-A-- MNDOAtom
// OrbitalProperties <--HAS-A-- LCGTO <--IS-A-- STO6G <--IS-A-- MNDO6G
public abstract class AtomFixed { // abstract class. Consists of an array of
	// LCGTO objects (ie. the valence orbitals), a coordinate array, an atomic
	// charge Z and a core charge Q.
	protected double[] coordinates;
	protected AtomProperties atomProperties;

	public AtomFixed(AtomProperties atomProperties, double[] coordinates) {
		// HAS AtomProperties and coordinates
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
