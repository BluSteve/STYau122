package runcycle.output;

import runcycle.input.RawMolecule;

import java.util.Arrays;

public class MoleculeOutput {
	public RawMolecule rawMolecule;
	public long time;
	public boolean isExpAvail;
	public double hf, dipole, ie, geomGradient, totalError;
	public ParamGradientOutput gradient;
	public double[][] hessian;

	@Override
	public String toString() {
		return "MoleculeOutput{" +
				"rawMolecule=" + rawMolecule +
				", time=" + time +
				", isExpAvail=" + isExpAvail +
				", hf=" + hf +
				", dipole=" + dipole +
				", ie=" + ie +
				", geomGradient=" + geomGradient +
				", totalError=" + totalError +
				", gradient=" + gradient +
				", hessian=" + Arrays.toString(hessian) +
				'}';
	}
}
