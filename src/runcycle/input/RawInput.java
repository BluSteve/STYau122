package runcycle.input;

import java.util.Arrays;

public class RawInput {
	public String trainingSet;
	public int[] atomTypes;
	public RawParams params;
	public RawMolecule[] molecules;

	@Override
	public String toString() {
		return "RawInput{" +
				"trainingSet='" + trainingSet + '\'' +
				", atomTypes=" + Arrays.toString(atomTypes) +
				", params=" + params +
				", molecules=" + Arrays.toString(molecules) +
				'}';
	}
}
