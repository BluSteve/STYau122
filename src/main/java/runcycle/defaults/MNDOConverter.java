package runcycle.defaults;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.defaults.MNDOAtom;
import nddo.structs.AtomProperties;
import runcycle.INDDOConverter;
import runcycle.structs.Atom;

public class MNDOConverter implements INDDOConverter {
	@Override
	public NDDOAtom convert(Atom atom, NDDOParams[] npMap) {
		return new MNDOAtom(AtomProperties.getAtoms()[atom.Z], atom.coords, npMap[atom.Z]);
	}
}
