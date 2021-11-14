package runcycle.input;

import java.util.Arrays;

public class RawParams {
	public double[][] nddoParams;
	public double[] lastHessian;
	public double[] lastGradient; // total gradients
	public double[] lastDir;

	@Override
	public String toString() {
		return "RawParams{" +
				"nddoParams=" + Arrays.toString(nddoParams) +
				", lastHessian=" + Arrays.toString(lastHessian) +
				", lastGradient=" + Arrays.toString(lastGradient) +
				", lastSearchDir=" + Arrays.toString(lastDir) +
				'}';
	}
}
