package runcycle.input;

public class RawParams {
	public final double[][] nddoParams;
	public final double[] lastDir;
	public final double[] lastGradient; // total gradients
	public final double[] lastHessian;

	public RawParams(double[][] nddoParams, double[] lastDir, double[] lastGradient, double[] lastHessian) {
		this.nddoParams = nddoParams;
		this.lastDir = lastDir;
		this.lastGradient = lastGradient;
		this.lastHessian = lastHessian;
	}
}
