package runcycle.output;

import runcycle.input.RawParams;

public class AggregateOutput {
	public double ttError;
	public long ttTime;
	public String inputHash, outputHash;
	public RawParams params;
	public MoleculeOutput[] mos;
}
