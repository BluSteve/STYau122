package runcycle.defaults;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.mndo.MNDOAtom;
import nddo.structs.AtomProperties;
import runcycle.NDDOConverter;
import runcycle.structs.Atom;

public class MNDOConverter implements NDDOConverter {
	private final NDDOParams[] npMap;

	public MNDOConverter(NDDOParams[] npMap) {
		this.npMap = npMap;
	}

	@Override
	public NDDOAtom convert(Atom atom) {
		return new MNDOAtom(AtomProperties.getAtoms()[atom.Z], atom.coords, npMap[atom.Z]);
	}

	@Override
	public NDDOAtom[] convert(Atom[] atoms) {
		NDDOAtom[] nas = new NDDOAtom[atoms.length];

		for (int i = 0; i < atoms.length; i++) {
			nas[i] = convert(atoms[i]);
		}

		return nas;
	}
}
