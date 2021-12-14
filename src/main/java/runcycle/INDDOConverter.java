package runcycle;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import runcycle.structs.Atom;

public interface INDDOConverter {
	NDDOAtom convert(Atom atom, NDDOParams[] npMap);
	NDDOAtom[] convert(Atom[] atoms, NDDOParams[] npMap);
}
