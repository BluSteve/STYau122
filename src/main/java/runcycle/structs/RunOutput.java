package runcycle.structs;

import runcycle.IMoleculeResult;

public final class RunOutput {
	public final String hash, inputHash, nextInputHash;
	public final long timeTaken;
	public final double ttError;
	public final double[] ttGradient;
	public final double[][] ttHessian;
	public final IMoleculeResult[] results;

	public final transient RunInput input, nextInput;

	public RunOutput(IMoleculeResult[] results, long timeTaken, double ttError, double[] ttGradient,
					 double[][] ttHessian, RunInput input, RunInput nextInput) {
		this.results = results;
		this.timeTaken = timeTaken;
		this.ttError = ttError;
		this.ttGradient = ttGradient;
		this.ttHessian = ttHessian;
		this.input = input;
		this.nextInput = nextInput;

		inputHash = input.hash;
		nextInputHash = nextInput.hash;

		hash = Serializer.getHash(this);
	}
}
