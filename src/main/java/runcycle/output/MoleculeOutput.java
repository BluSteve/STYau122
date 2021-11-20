package runcycle.output;

import runcycle.input.RawMolecule;

public class MoleculeOutput {
	public RawMolecule rawMolecule;
	public long time;
	public double hf, dipole, ie, geomGradient, totalError;
	public ParamGradientOutput gradient;
	public double[][] hessian;
}
