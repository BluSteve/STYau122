package runcycle.structs;

public class LastRunInfo {
	public double error, trustRadius, stepSize, expectedChange;

	@Override
	public String toString() {
		return "LastRunInfo{" +
				"error=" + error +
				", trustRadius=" + trustRadius +
				", stepSize=" + stepSize +
				", expectedChange=" + expectedChange +
				'}';
	}
}
