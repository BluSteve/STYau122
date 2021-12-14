package frontend;

import runcycle.structs.Params;
import runcycle.structs.RunnableMolecule;

public class RawInput {
	public String hash, trainingSet;
	public int nMolecules;
	public int[] atomTypes;
	public int[][] neededParams;
	public Params params;
	public RunnableMolecule[] molecules;
}
