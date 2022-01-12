package runcycle.structs;

import nddo.Constants;
import nddo.NDDOParams;

import java.util.ArrayList;
import java.util.List;

public class InputInfo {
	public final int[] atomTypes;
	public final int[][] neededParams; // todo make this a map
	public final NDDOParams[] npMap;

	public InputInfo(int[] atomTypes, int[][] neededParams, int[] paramsAtomTypes, double[][] params) {
		this.atomTypes = atomTypes;
		this.neededParams = neededParams;
		this.npMap = constructNpMap(paramsAtomTypes, params);
	}

	public InputInfo(int[] atomTypes, int[][] neededParams, double[][] paramsMap) {
		this.atomTypes = atomTypes;
		this.neededParams = neededParams;

		npMap = new NDDOParams[paramsMap.length];

		for (int i = 0; i < paramsMap.length; i++) {
			if (paramsMap[i].length > 0) npMap[i] = new NDDOParams(paramsMap[i]);
		}
	}

	public InputInfo(int[] atomTypes, int[][] neededParams, NDDOParams[] npMap) {
		this.atomTypes = atomTypes;
		this.neededParams = neededParams;
		this.npMap = npMap;
	}

	private static NDDOParams[] constructNpMap(int[] atomTypes, double[][] params) {
		if (atomTypes.length > params.length) throw new NullPointerException("Not enough params for atomTypes!");

		NDDOParams[] result = new NDDOParams[Constants.MAX_ATOM_NUM];

		for (int i = 0; i < atomTypes.length; i++) {
			result[atomTypes[i]] = new NDDOParams(params[i]);
		}

		return result;
	}

	public double[][] getParams() {
		List<double[]> res = new ArrayList<>();

		for (int i = 0; i < Constants.MAX_ATOM_NUM; i++) {
			if (npMap[i] != null) res.add(npMap[i].params.clone());
		}

		return res.toArray(new double[0][]);
	}

	public double[][] getParamsMap() {
		double[][] res = new double[Constants.MAX_ATOM_NUM][];

		for (int i = 0; i < Constants.MAX_ATOM_NUM; i++) {
			res[i] = npMap[i] != null ? npMap[i].params.clone(): new double[0];
		}

		return res;
	}
}
