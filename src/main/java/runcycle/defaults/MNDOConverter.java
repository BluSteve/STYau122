package runcycle.defaults;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.mndo.MNDOAtom;
import nddo.structs.AtomProperties;
import runcycle.INDDOConverter;
import runcycle.structs.Atom;

public class MNDOConverter implements INDDOConverter {
	@Override
	public NDDOAtom convert(Atom atom, NDDOParams[] npMap) {
		return new MNDOAtom(AtomProperties.getAtoms()[atom.Z], atom.coords, npMap[atom.Z]);
	}

	@Override
	public NDDOAtom[] convert(Atom[] atoms, NDDOParams[] npMap) {
		NDDOAtom[] nas = new NDDOAtom[atoms.length];

		for (int i = 0; i < atoms.length; i++) {
			nas[i] = convert(atoms[i], npMap);
		}

		return nas;
	}
}
