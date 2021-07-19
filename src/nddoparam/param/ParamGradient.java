package nddoparam.param;

import nddoparam.Solution;
import org.jblas.DoubleMatrix;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;

public abstract class ParamGradient {
	protected static final double LAMBDA = 1E-7;
	protected Solution s, sPrime, sExp;
	protected ParamErrorFunction e;
	protected boolean isExpAvail, analytical;
	protected double[] datum;
	protected double[][] HFDerivs, dipoleDerivs, IEDerivs, geomDerivs,
			totalGradients;
	protected DoubleMatrix[][] densityDerivs, xLimited, xComplementary, xForIE,
			coeffDerivs, responseDerivs, fockDerivs;
	protected DoubleMatrix[][][] staticDerivs;

	public ParamGradient(Solution s, double[] datum, Solution sExp,
						 boolean analytical) {
		this.s = s;
		this.datum = datum;
		this.sExp = sExp;
		this.analytical = analytical;
		initializeArrays();
	}

	public static double[][] depad(double[][] derivs, int[][] neededParams,
								   int[] uniqueZs) {
		double[][] res = new double[derivs.length][0];
		for (int i = 0; i < derivs.length; i++) {
			double[] depadded =
					new double[neededParams[uniqueZs[i]].length];
			int u = 0;
			for (int j = 0; j < derivs[i].length; j++) {
				for (int p : neededParams[uniqueZs[i]]) {
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

	// also pads atomTypes that this molecule doesn't contain with zeros, for
	// consistency's sake.
	public static double[] combine(double[][] derivs, int[] atomTypes,
								   int[][] neededParams, int[] uniqueZs,
								   int[][] moleculeNP, boolean isDepad) {
		if (isDepad) derivs = depad(derivs, moleculeNP, uniqueZs);

		// confusingly, this is a different kind of padding from the above
		// depad method. depad removes the zeros which represent
		// non-differentiated parameters, but this adds padding for the atoms
		// that are in the training set but not this particular molecule.
		double[][] paddedDerivs = new double[atomTypes.length][];
		for (int i = 0; i < atomTypes.length; i++) {
			for (int j = 0; j < uniqueZs.length; j++) {
				if (atomTypes[i] == uniqueZs[j]) {
					paddedDerivs[i] = derivs[j];
				}
				else {
					paddedDerivs[i] =
							new double[neededParams[i].length];
				}
			}
		}

		ArrayList<Double> a = new ArrayList<>();
		for (double[] deriv : paddedDerivs) {
			for (double d : deriv) {
				a.add(d);
			}
		}
		double[] res = new double[a.size()];
		for (int i = 0; i < res.length; i++) res[i] = a.get(i);

		return res;
	}

	public void compute() {
		totalGradients =
				new double[s.getUniqueZs().length][Solution.maxParamNum];
		if (analytical && (datum[1] != 0 || datum[2] != 0))
			computeBatchedDerivs(0, 0);

		// if geom finite difference computations are not needed, overhead
		// for multithreading exceeds potential gain.
		if (!analytical || isExpAvail) {
			List<RecursiveAction> subtasks = new ArrayList<>();

			for (int Z = 0; Z < s.getUniqueZs().length; Z++) {
				for (int paramNum : s.getNeededParams()[s.getUniqueZs()[Z]]) {
					int finalZ = Z;
					subtasks.add(new RecursiveAction() {
						@Override
						protected void compute() {
							computeGradient(finalZ, paramNum);
						}
					});
				}
			}
			ForkJoinTask.invokeAll(subtasks);
		}
		else {
			for (int Z = 0; Z < s.getUniqueZs().length; Z++) {
				for (int paramNum : s.getNeededParams()[s.getUniqueZs()[Z]]) {
					computeGradient(Z, paramNum);
				}
			}
		}
	}

	protected void computeGradient(int Z, int paramNum) {
		Solution sPrime = null;
		if (!analytical) {
			sPrime = constructSPrime(Z, paramNum);
		}

		computeHFDeriv(Z, paramNum, sPrime);
		if (datum[1] != 0 && datum[2] != 0) {
			computeDipoleDeriv(Z, paramNum, true, sPrime);
			computeIEDeriv(Z, paramNum, sPrime);
		}
		else if (datum[1] != 0) {
			computeDipoleDeriv(Z, paramNum, true, sPrime);
		}
		else if (datum[2] != 0) {
			computeDipoleDeriv(Z, paramNum, false, sPrime);
			computeIEDeriv(Z, paramNum, sPrime);
		}

		if (isExpAvail) computeGeomDeriv(Z, paramNum);
	}

	protected void computeGeomDeriv(int Z, int paramNum) {
		Solution sExpPrime = constructSExpPrime(Z, paramNum);
		double sum = 0;
		double d;
		for (int i = 0; i < sExpPrime.atoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				d = findGrad(sExpPrime, i, j);
				sum += d * d;
			}
		}
		double geomGradient = 627.5 * Math.sqrt(sum);
		geomDerivs[Z][paramNum] = 1 / LAMBDA * (geomGradient - e.geomGradient);
		totalGradients[Z][paramNum] +=
				0.000098 * e.geomGradient * geomDerivs[Z][paramNum];
	}

	protected void errorFunctionRoutine() {
		if (datum[1] != 0) e.addDipoleError(datum[1]);
		if (datum[2] != 0) e.addIEError(datum[2]);

		if (this.sExp != null) {
			isExpAvail = true;
			geomDerivs = new double[s.getUniqueZs().length]
					[Solution.maxParamNum];
			e.createExpGeom(this.sExp);
			e.addGeomError();
		}
	}

	protected void initializeArrays() {
		int atomLength = s.getUniqueZs().length;
		int paramLength = Solution.maxParamNum;
		totalGradients = new double[atomLength][paramLength];
		HFDerivs = new double[atomLength][paramLength];
		if (datum[1] != 0) dipoleDerivs = new double[atomLength][paramLength];
		if (datum[2] != 0) IEDerivs = new double[atomLength][paramLength];

		if (analytical && (datum[1] != 0 || datum[2] != 0)) {
			densityDerivs = new DoubleMatrix[atomLength][paramLength];
			staticDerivs = new DoubleMatrix[atomLength][2][paramLength];
			xLimited = new DoubleMatrix[atomLength][paramLength];
			if (datum[2] != 0) {
				xComplementary = new DoubleMatrix[atomLength][paramLength];
				xForIE = new DoubleMatrix[atomLength][paramLength];
				coeffDerivs = new DoubleMatrix[atomLength][paramLength];
				responseDerivs = new DoubleMatrix[atomLength][paramLength];
				fockDerivs = new DoubleMatrix[atomLength][paramLength];
			}
		}

	}

	protected abstract Solution constructSPrime(int Z, int paramNum);

	protected abstract void computeBatchedDerivs(int firstZIndex,
												 int firstParamIndex);

	protected abstract void computeHFDeriv(int Z, int paramNum,
										   Solution sPrime);

	protected abstract void computeDipoleDeriv(int Z, int paramNum,
											   boolean full,
											   Solution sPrime);

	protected abstract void computeIEDeriv(int Z, int paramNum,
										   Solution sPrime);

	protected abstract Solution constructSExpPrime(int Z, int paramNum);

	protected abstract double findGrad(Solution sExpPrime, int i, int j);

	public ParamErrorFunction getE() {
		return this.e;
	}

	public Solution getS() {
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

	@Deprecated
	public double getAnalyticalError() {
		this.compute();
		double[][] a = new double[totalGradients.length][0];
		for (int i = 0; i < totalGradients.length; i++)
			a[i] = totalGradients[i].clone();
		analytical = !analytical;
		this.compute();
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
