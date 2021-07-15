package runcycle.input;

import java.util.Arrays;

public class RawInput {
	public String model, trainingSet;
	public int nMolecules;
	public int[] atomTypes;
	public RawParams params;
	public RawMolecule[] molecules;

	@Override
	public String toString() {
		return "RawInput{" +
				"trainingSet='" + trainingSet + '\'' +
				", nMolecules=" + nMolecules +
				", atomTypes=" + Arrays.toString(atomTypes) +
				", params=" + params +
				", molecules=" + Arrays.toString(molecules) +
				'}';
	}
}
