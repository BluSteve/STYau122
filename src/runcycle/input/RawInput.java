package runcycle.input;

import java.util.Arrays;

public class RawInput {
	public String hash, model, trainingSet;
	public int nMolecules;
	public int[] atomTypes;
	public int[][] neededParams;
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
