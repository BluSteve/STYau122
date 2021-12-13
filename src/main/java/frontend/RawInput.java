package frontend;

import nddo.structs.Model;
import runcycle.structs.RunnableMolecule;

public class RawInput {
	public String hash, trainingSet;
	public Model model; // todo don't hardcode model
	public int nMolecules;
	public int[] atomTypes;
	public int[][] neededParams;
	public RawParams params;
	public RunnableMolecule[] molecules;
}
