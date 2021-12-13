package runcycle.input;

public class InputInfo {
	public final int[] atomTypes;
	public final int[][] neededParams;
	public final RawParams params;

	public InputInfo(int[] atomTypes, int[][] neededParams, RawParams params) {
		this.atomTypes = atomTypes;
		this.neededParams = neededParams;
		this.params = params;
	}
}
