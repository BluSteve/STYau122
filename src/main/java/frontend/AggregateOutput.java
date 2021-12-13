package frontend;

import runcycle.output.MoleculeOutput;

public class AggregateOutput {
	public double ttError;
	public long ttTime;
	public String inputHash, outputHash;
	public RawParams params;
	public MoleculeOutput[] mos;
}
