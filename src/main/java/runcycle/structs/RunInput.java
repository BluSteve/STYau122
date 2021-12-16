package runcycle.structs;

public class RunInput {
	public final InputInfo info; // info can vary independently of molecules
	public final RunnableMolecule[] molecules;

	public RunInput(InputInfo info, RunnableMolecule[] molecules) {
		this.info = info;
		this.molecules = molecules;
	}
}
