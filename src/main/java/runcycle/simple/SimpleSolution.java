package runcycle.simple;

import nddo.NDDOAtom;
import nddo.solution.Solution;
import runcycle.State;
import runcycle.structs.Atom;
import runcycle.structs.InputInfo;
import runcycle.structs.RunnableMolecule;

public class SimpleSolution {
	public Solution of(InputInfo info, Atom[] atoms, boolean restricted, int charge, int mult) {
		NDDOAtom[] nddoAtoms = State.getConverter().convert(atoms, info.npMap);

		RunnableMolecule.RMBuilder builder = new RunnableMolecule.RMBuilder();
		builder.atoms = atoms;
		builder.restricted = restricted;
		builder.charge = charge;
		builder.mult = mult;

		RunnableMolecule rm = builder.build(info.atomTypes, info.neededParams, info.npMap);

		return Solution.of(rm, nddoAtoms);
	}
}
