package runcycle.structs;

public final class RunInput {
	public String hash;
	public InputInfo info; // info can vary independently of molecules
	public RunnableMolecule[] molecules;

	public RunInput(InputInfo info, RunnableMolecule[] molecules) {
		this.info = info;
		this.molecules = molecules;

		hash = Serializer.getHash(this);
	}
}
