package runcycle;

import nddo.NDDOAtom;
import runcycle.structs.Atom;

public interface NDDOConverter {
	NDDOAtom convert(Atom atom);
	NDDOAtom[] convert(Atom[] atoms);
}
