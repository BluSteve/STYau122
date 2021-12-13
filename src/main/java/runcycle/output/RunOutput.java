package runcycle.output;

import runcycle.input.InputInfo;

public class RunOutput {
	public final InputInfo nextRunInfo;
	public final MoleculeOutput[] mos;

	public RunOutput(InputInfo nextRunInfo, MoleculeOutput[] mos) {
		this.nextRunInfo = nextRunInfo;
		this.mos = mos;
	}
}