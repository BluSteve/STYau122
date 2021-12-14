package runcycle.structs;

import nddo.NDDOParams;

public class Params {
	public final NDDOParams[] npMap;
	public final double[] lastDir;
	public final double[] lastGradient; // total gradients
	public final double[] lastHessian;

	public Params(NDDOParams[] npMap, double[] lastDir, double[] lastGradient, double[] lastHessian) {
		this.npMap = npMap;
		this.lastDir = lastDir;
		this.lastGradient = lastGradient;
		this.lastHessian = lastHessian;
	}

	/**
	 *
	 * @param totalParamLength Combined length of all trainable parameters.
	 */
	public Params(NDDOParams[] npMap, int totalParamLength) {
		this.npMap = npMap;
		this.lastDir = new double[totalParamLength];
		this.lastGradient = new double[totalParamLength];
		this.lastHessian = new double[(totalParamLength + 1) * totalParamLength / 2];
	}
}
