package runcycle;

import nddo.NDDOAtom;
import runcycle.structs.Atom;

public interface INDDOConverter {
	NDDOAtom convert(Atom atom);
	NDDOAtom[] convert(Atom[] atoms);
}
