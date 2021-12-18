package runcycle.structs;

import runcycle.IMoleculeResult;

public final class RunOutput { // non-transient fields must be for reference purposes only
	public final String hash;
	public final InputInfo nextRunInfo;
	public final IMoleculeResult[] results;
	public final long timeTaken;
	private final transient RunInput input; // what made it
	private transient RunInput nextInput; // what it made

	public RunOutput(RunInput input, InputInfo nextRunInfo, IMoleculeResult[] results, long timeTaken) {
		this.input = input;
		this.nextRunInfo = nextRunInfo;
		this.results = results;
		this.timeTaken = timeTaken;

		hash = Serializer.getHash(this);
	}

	public RunInput getInput() {
		return input;
	}

	public RunInput getNextInput() {
		if (nextInput == null) {
			RunnableMolecule[] nextRunRms = new RunnableMolecule[results.length];
			for (int i = 0; i < nextRunRms.length; i++) {
				nextRunRms[i] = results[i].getUpdatedRm();
			}

			nextInput = new RunInput(nextRunInfo, nextRunRms);
		}

		return nextInput;
	}
}