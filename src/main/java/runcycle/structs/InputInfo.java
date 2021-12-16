package runcycle.structs;

import nddo.Constants;
import nddo.NDDOParams;

public class InputInfo {
	public final int[] atomTypes;
	public final int[][] neededParams;
	public final NDDOParams[] npMap;

	public InputInfo(int[] atomTypes, int[][] neededParams, double[][] params) {
		this.atomTypes = atomTypes;
		this.neededParams = neededParams;
		this.npMap = constructNpMap(atomTypes, params);
	}

	public InputInfo(int[] atomTypes, int[][] neededParams, NDDOParams[] npMap) {
		this.atomTypes = atomTypes;
		this.neededParams = neededParams;
		this.npMap = npMap;
	}

	private static NDDOParams[] constructNpMap(int[] Zs, double[][] params) {
		NDDOParams[] result = new NDDOParams[Constants.maxAtomNum];

		for (int i = 0; i < Zs.length; i++) {
			result[Zs[i]] = new NDDOParams(params[i]);
		}

		return result;
	}

	public double[][] getParams() {
		double[][] res = new double[atomTypes.length][];

		for (int i = 0; i < atomTypes.length; i++) {
			int Z = atomTypes[i];
			res[i] = npMap[Z].toArray();
		}

		return res;
	}
}
