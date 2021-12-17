package runcycle.structs;

import runcycle.Serializer;

public final class RunInput {
	public final String hash;
	public final InputInfo info; // info can vary independently of molecules
	public final RunnableMolecule[] molecules;

	public RunInput(InputInfo info, RunnableMolecule[] molecules) {
		this.info = info;
		this.molecules = molecules;

		hash = Serializer.getHash(this);
	}
}
