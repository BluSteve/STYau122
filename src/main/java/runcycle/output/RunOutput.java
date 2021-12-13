package runcycle.output;

import runcycle.input.InputInfo;

public interface RunOutput {
	InputInfo getNextRunInfo();
	MoleculeOutput[] getMOs();
}
