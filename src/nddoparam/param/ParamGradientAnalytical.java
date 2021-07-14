package nddoparam.param;

import nddoparam.NDDOSolution;
import org.jblas.DoubleMatrix;

public abstract class ParamGradientAnalytical implements ErrorGettable {
	protected static final double LAMBDA = 1E-7;
	protected NDDOSolution s, sPrime, sExpPrime, sExp;
	protected ParamErrorFunction e;
	protected String kind;
	protected boolean isExpAvail, analytical;
	protected double[] datum;
	protected double[][] HFDerivs, dipoleDerivs, IEDerivs, geomDerivs, totalGradients;
	protected DoubleMatrix[][] densityDerivs, xLimited, xComplementary, xForIE,
			coeffDerivs, responseDerivs, fockDerivs;
	protected DoubleMatrix[][][] staticDerivs;

	public ParamGradientAnalytical(NDDOSolution s, String kind, double[] datum,
								   NDDOSolution sExp, boolean analytical) {
		this.s = s;
		this.kind = kind;
		this.datum = datum;
		this.sExp = sExp;
		this.analytical = analytical;
	}

	protected void errorFunctionRoutine() {
		if (datum[1] != 0) e.addDipoleError(datum[1]);
		if (datum[2] != 0) e.addIEError(datum[2]);

		if (this.sExp != null) {
			isExpAvail = true;
			geomDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
			e.createExpGeom(this.sExp);
			e.addGeomError();
		}
	}

	protected void initializeArrays() {
		totalGradients = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
		switch (kind) {
			case "a":
				HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
				break;
			case "b":
				HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
				dipoleDerivs =
						new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

				if (this.analytical) {
					densityDerivs = new DoubleMatrix[s
							.getUniqueZs().length][NDDOSolution.maxParamNum];
					staticDerivs = new DoubleMatrix[s
							.getUniqueZs().length][NDDOSolution.maxParamNum][2];
					xLimited = new DoubleMatrix[s
							.getUniqueZs().length][NDDOSolution.maxParamNum];
				}
				break;
			case "c":
				HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
				IEDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

				if (analytical) initializeIntermediates();
				break;
			case "d":
				HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
				dipoleDerivs =
						new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
				IEDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

				if (analytical) initializeIntermediates();
				break;
		}
	}

	private void initializeIntermediates() {
		densityDerivs =
				new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
		staticDerivs =
				new DoubleMatrix[s.getUniqueZs().length][2][NDDOSolution.maxParamNum];
		xLimited = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];

		xComplementary =
				new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
		xForIE = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
		coeffDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
		responseDerivs =
				new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
		fockDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
	}

	public void computeGradients() {
		totalGradients = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
		if (analytical && (kind.equals("b") || kind.equals("c") || kind.equals("d")))
			computeBatchedDerivs(0, 0);
		for (int Z = 0; Z < s.getUniqueZs().length; Z++) {
			for (int paramNum : s.getNeededParams()[s.getUniqueZs()[Z]]) {
				computeGradient(Z, paramNum);
			}
		}
	}

	public void computeGradient(int Z, int paramNum) {
		if (!analytical) constructSPrime(Z, paramNum);
		switch (kind) {
			case "a":
				computeHFDeriv(Z, paramNum);
				break;
			case "b":
				computeDipoleDeriv(Z, paramNum, true);
				break;
			case "c":
				computeDipoleDeriv(Z, paramNum, false);
				computeIEDeriv(Z, paramNum);
				break;
			case "d":
				computeDipoleDeriv(Z, paramNum, true);
				computeIEDeriv(Z, paramNum);
				break;
		}
		if (isExpAvail) computeGeomDeriv(Z, paramNum);
	}

	protected void computeGeomDeriv(int Z, int paramNum) {
		constructSExpPrime(Z, paramNum);
		double sum = 0;
		double d;
		for (int i = 0; i < sExpPrime.atoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				d = findGrad(i, j);
				sum += d * d;
			}
		}
		double geomGradient = 627.5 * Math.sqrt(sum);
		geomDerivs[Z][paramNum] = 1 / LAMBDA * (geomGradient - e.geomGradient);
		totalGradients[Z][paramNum] +=
				0.000098 * e.geomGradient * geomDerivs[Z][paramNum];
	}

	protected abstract void constructSPrime(int Z, int paramNum);

	protected abstract void computeBatchedDerivs(int firstZIndex, int firstParamIndex);

	protected abstract void computeHFDeriv(int Z, int paramNum);

	protected abstract void computeDipoleDeriv(int Z, int paramNum, boolean full);

	protected abstract void computeIEDeriv(int Z, int paramNum);

	protected abstract void constructSExpPrime(int Z, int paramNum);

	protected abstract double findGrad(int i, int j);

	public ParamErrorFunction getE() {
		return this.e;
	}

	public NDDOSolution getS() {
		return s;
	}

	public double[][] getHFDerivs() {
		return HFDerivs;
	}

	public double[][] getDipoleDerivs() {
		return dipoleDerivs;
	}

	public double[][] getIEDerivs() {
		return IEDerivs;
	}

	public double[][] getGeomDerivs() {
		return geomDerivs;
	}

	public double[][] getTotalGradients() {
		return totalGradients;
	}

	public boolean isAnalytical() {
		return analytical;
	}

	public void setAnalytical(boolean analytical) {
		this.analytical = analytical;
	}

	public double[][] depad(double[][] derivs) {
		double[][] res = new double[derivs.length][0];
		for (int i = 0; i < derivs.length; i++) {
			double[] depadded =
					new double[s.getNeededParams()[s.getUniqueZs()[i]].length];
			int u = 0;
			for (int j = 0; j < derivs[i].length; j++) {
				for (int p : s.getNeededParams()[s.getUniqueZs()[i]]) {
					if (j == p) {
						depadded[u] = derivs[i][j];
						u++;
						break;
					}
				}
			}
			res[i] = depadded;
		}
		return res;
	}

	public double getAnalyticalError() {
		this.computeGradients();
		double[][] a = new double[totalGradients.length][0];
		for (int i = 0; i < totalGradients.length; i++) a[i] = totalGradients[i].clone();
		analytical = !analytical;
		this.computeGradients();
		double sum = 0;
		for (int i = 0; i < totalGradients.length; i++) {
			for (int j = 0; j < totalGradients[0].length; j++) {
				sum += (totalGradients[i][j] - a[i][j]) *
						(totalGradients[i][j] - a[i][j]);
			}
		}
		analytical = !analytical;
		totalGradients = a;
		return sum;
	}
}
