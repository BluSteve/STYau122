package runcycle.structs;

import runcycle.IMoleculeResult;

public class RunOutput {
	public final InputInfo nextRunInfo;
	public final IMoleculeResult[] results;
	public final long timeTaken;

	public RunOutput(InputInfo nextRunInfo, IMoleculeResult[] results, long timeTaken) {
		this.nextRunInfo = nextRunInfo;
		this.results = results;
		this.timeTaken = timeTaken;
	}
}