package runcycle.input;

import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import scf.AtomHandler;

import java.util.Arrays;

public class RawMolecule {
	public int index;
	public String name;
	public boolean isUsing, restricted;
	public int charge, mult, nElectrons;
	public double[] datum;
	public int[] uniqueZs;
	public RawAtom[] atoms, expGeom;

	public static MNDOAtom[] toMNDOAtoms(RawAtom[] atoms,
										 MNDOParams[] mndoParams) {
		MNDOAtom[] res = new MNDOAtom[atoms.length];
		for (int i = 0; i < atoms.length; i++) {
			res[i] = atoms[i].toMNDOAtom(
					mndoParams[AtomHandler.atoms[atoms[i].Z].getIndex()]);
		}
		return res;
	}

	@Override
	public String toString() {
		return "RawMolecule{" +
				"index=" + index +
				", name='" + name + '\'' +
				", restricted=" + restricted +
				", charge=" + charge +
				", mult=" + mult +
				", nElectrons=" + nElectrons +
				", datum=" + Arrays.toString(datum) +
				", uniqueZs=" + Arrays.toString(uniqueZs) +
				", atoms=" + Arrays.toString(atoms) +
				", expGeom=" + Arrays.toString(expGeom) +
				'}';
	}
}
