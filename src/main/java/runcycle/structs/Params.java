package runcycle.structs;

import nddo.Constants;
import nddo.NDDOParams;

public class Params {
	public final NDDOParams[] npMap;

	public Params(int[] Zs, double[][] params) {
		this.npMap = constructNpMap(Zs, params);
	}

	public Params(NDDOParams[] npMap) {
		this.npMap = npMap;
	}

	private static NDDOParams[] constructNpMap(int[] Zs, double[][] params) {
		NDDOParams[] result = new NDDOParams[Constants.maxAtomNum];

		for (int i = 0; i < Zs.length; i++) {
			result[Zs[i]] = new NDDOParams(params[i]);
		}

		return result;
	}
}
