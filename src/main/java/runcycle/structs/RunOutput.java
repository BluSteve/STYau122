package runcycle.structs;

import runcycle.IMoleculeResult;

public final class RunOutput {
	public String hash, inputHash, nextInputHash;
	public long timeTaken;
	public double ttError;
	public double finalLambda;
	public double[] ttGradient;
	public double[][] ttHessian;
	public IMoleculeResult[] results;
	public transient RunInput input, nextInput;

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
