package examples.am1;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.structs.AtomProperties;
import runcycle.INDDOConverter;
import runcycle.structs.Atom;

public class AM1Converter implements INDDOConverter {
	@Override
	public NDDOAtom convert(Atom atom, NDDOParams[] npMap) {
		return new AM1Atom(AtomProperties.getAtoms()[atom.Z], atom.coords.clone(), npMap[atom.Z]);
	}
}
