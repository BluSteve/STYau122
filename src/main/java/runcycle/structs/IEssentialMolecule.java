package runcycle.structs;

public interface IEssentialMolecule {
	int getCharge();

	int getMult();

	boolean isRestricted();

	Atom[] getAtoms();

	Atom[] getExpGeom();

	double[] getDatum();

	int[] getMats();

	int[][] getMnps();

	String debugName();

	String getName();

	int getIndex();
}
