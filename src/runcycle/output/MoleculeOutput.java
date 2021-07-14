package runcycle.output;

public class MoleculeOutput {
	public int index;
	public String name;
	public long time;
	public double[] datum;
	public double hf, dipole, ie, geomGradient, totalError;
	public ParamGradientOutput gradient;
	public double[][] hessian;
}
