package runcycle;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import runcycle.structs.Atom;

public interface INDDOConverter {
	NDDOAtom convert(Atom atom, NDDOParams[] npMap);

	default NDDOAtom[] convert(Atom[] atoms, NDDOParams[] npMap) {
		NDDOAtom[] nas = new NDDOAtom[atoms.length];

		for (int i = 0; i < atoms.length; i++) {
			nas[i] = convert(atoms[i], npMap);
		}

		return nas;
	}
}
