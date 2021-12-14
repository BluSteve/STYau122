package runcycle.structs;

public class InputInfo {
	public final int[] atomTypes;
	public final int[][] neededParams;
	public final Params params;

	public InputInfo(int[] atomTypes, int[][] neededParams, Params params) {
		this.atomTypes = atomTypes;
		this.neededParams = neededParams;
		this.params = params;
	}
}
