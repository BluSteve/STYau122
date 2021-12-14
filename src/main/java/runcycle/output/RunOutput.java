package runcycle.output;

import runcycle.IMoleculeResult;
import runcycle.input.InputInfo;

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