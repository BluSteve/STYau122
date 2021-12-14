package runcycle.structs;

import nddo.Constants;
import nddo.NDDOParams;

public class Params {
	public final NDDOParams[] npMap;
	public final double[] lastDir;
	public final double[] lastGradient; // total gradients
	public final double[] lastHessian;

	public Params(int[] Zs, double[][] params, double[] lastDir, double[] lastGradient, double[] lastHessian) {
		this.npMap = constructNpMap(Zs, params);
		this.lastDir = lastDir;
		this.lastGradient = lastGradient;
		this.lastHessian = lastHessian;
	}

	/**
	 *
	 * @param totalParamLength Combined length of all trainable parameters.
	 */
	public Params(int[] Zs, double[][] params, int totalParamLength) {
		this.npMap = constructNpMap(Zs, params);
		this.lastDir = new double[totalParamLength];
		this.lastGradient = new double[totalParamLength];
		this.lastHessian = new double[(totalParamLength + 1) * totalParamLength / 2];
	}

	private static NDDOParams[] constructNpMap(int[] Zs, double[][] params) {
		NDDOParams[] result = new NDDOParams[Constants.maxAtomNum];

		for (int i = 0; i < Zs.length; i++) {
			result[Zs[i]] = new NDDOParams(params[i]);
		}

		return result;
	}
}
