package runcycle.output;

public class MoleculeOutput {
	int index;
	String name;
	long time;
	double[] datum;
	double hf, dipole, ie, geomGradient, total;
	ParamGradientOutput gradient;
	double[][] hessian;
}
