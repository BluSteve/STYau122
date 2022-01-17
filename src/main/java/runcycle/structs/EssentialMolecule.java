package runcycle.structs;

public class EssentialMolecule implements  IEssentialMolecule{
	public int index;
	public String name, debugName;
	public boolean restricted;
	public int charge, mult;
	public double[] datum;
	public int[] mats;
	public int[][] mnps;
	public Atom[] atoms, expGeom;

	public EssentialMolecule() {
	}

	public EssentialMolecule(IEssentialMolecule iem) {
		index = iem.getIndex();
		atoms = iem.getAtoms();
		expGeom = iem.getExpGeom();
		datum = iem.getDatum();
		name = iem.getName();
		mats = iem.getMats();
		mnps = iem.getMnps();
		charge = iem.getCharge();
		mult = iem.getMult();
		restricted = iem.isRestricted();
		debugName = debugName();
	}

	public boolean isRestricted() {
		return restricted;
	}

	@Override
	public int getCharge() {
		return charge;
	}

	@Override
	public int getMult() {
		return mult;
	}

	@Override
	public String debugName() {
		if (debugName == null)
			debugName = String.format("%03d-%s_%d_%s", index, name, charge, restricted ? "RHF" : "UHF");

		return debugName;
	}

	@Override
	public int getIndex() {
		return index;
	}

	@Override
	public Atom[] getAtoms() {
		return atoms;
	}

	@Override
	public Atom[] getExpGeom() {
		return expGeom;
	}

	@Override
	public double[] getDatum() {
		return datum;
	}

	@Override
	public String getName() {
		return name;
	}

	@Override
	public int[] getMats() {
		return mats;
	}

	@Override
	public int[][] getMnps() {
		return mnps;
	}
}
