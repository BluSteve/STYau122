package frontend.json;

import runcycle.structs.Params;

public class AggregateOutput {
	public double ttError;
	public long ttTime;
	public String inputHash, outputHash;
	public Params params;
	public MoleculeOutput[] mos;
}
