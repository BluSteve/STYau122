package runcycle.structs;

public final class RunInput {
	public String hash;
	public InputInfo info; // info can vary independently of molecules
	public RunnableMolecule[] molecules;
	public LastRunInfo lastRunInfo;

	public RunInput(InputInfo info, RunnableMolecule[] molecules, LastRunInfo lastRunInfo) {
		this.info = info;
		this.molecules = molecules;
		this.lastRunInfo = lastRunInfo;

		hash = Serializer.getHash(this);
	}

	public RunInput(InputInfo info, RunnableMolecule[] molecules) {
		this(info, molecules, null);
	}
}
