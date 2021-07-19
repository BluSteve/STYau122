package nddoparam.param;

import nddoparam.Solution;
import nddoparam.SolutionR;
import nddoparam.SolutionU;
import org.jblas.DoubleMatrix;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;

public abstract class ParamGradient {
	protected static final double LAMBDA = 1E-7;
	protected Solution s, sExp;
	protected ParamErrorFunction e;
	protected boolean isExpAvail, analytical;
	protected double[] datum;
	protected double[][] HFDerivs, dipoleDerivs, IEDerivs, geomDerivs,
			totalGradients;
	protected DoubleMatrix[][] densityDerivs, xLimited, xComplementary, xForIE,
			coeffDerivs, responseDerivs, fockDerivs;
	protected DoubleMatrix[][][] staticDerivs;

	protected ParamGradient(Solution s, double[] datum, Solution sExp,
							boolean analytical) {
		this.s = s;
		this.datum = datum;
		this.sExp = sExp;
		this.analytical = analytical;
		e = ParamErrorFunction.of(s, datum[0]);

		initializeArrays();
		initializeErrorFunction();
	}

	/**
	 * Get gradients from basic ingredients.
	 *
	 * @param s     Primary Solution object.
	 * @param datum Reference data of Hf, dipole, I.E.
	 * @param sExp  Optional experimental geometry Solution object, may be
	 *                 null.
	 * @return Uncomputed ParamGradient object, either restricted or not
	 * depending on type of s and sExp.
	 */
	public static ParamGradient of(Solution s, double[] datum, Solution sExp) {
		if (s instanceof SolutionR)
			return new ParamGradientR((SolutionR) s, datum, (SolutionR) sExp,
					true);
		assert s instanceof SolutionU;
		return new ParamGradientU((SolutionU) s, datum, (SolutionU) sExp,
				false);
	}

	/**
	 * Depads any 2d deriv array based on arrays of atom types and their
	 * corresponding params required.
	 *
	 * @param derivs       2d deriv array. (e.g. HFDerivs, totalGradients, ...)
	 * @param neededParams Array of params needed in the same order as atom
	 *                     types. Assumes neededParams.length = derivs.length.
	 * @return Depadded 2d-array of the input derivs.
	 */
	public static double[][] depad(double[][] derivs, int[][] neededParams) {
		assert derivs.length == neededParams.length;
		double[][] res = new double[derivs.length][0];
		for (int i = 0; i < derivs.length; i++) {
			double[] depadded = new double[neededParams[i].length];
			int u = 0;
			for (int j = 0; j < derivs[i].length; j++) {
				for (int p : neededParams[i]) {
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

	/**
	 * Depads, pads, and flattens any 2d deriv array to make for easier
	 * comparison when outputted.
	 * <p/>
	 * Confusingly, padding is very different from depadding. depad removes
	 * the zeros which represent non-differentiated parameters, but this adds
	 * padding for the atoms that are in the training set but not this
	 * particular molecule.
	 *
	 * @param derivs       2d deriv array. (e.g. HFDerivs, totalGradients, ...)
	 * @param atomTypes    Array of atom types.
	 * @param neededParams Array of params needed in the same order as atom
	 *                     types.
	 * @param moleculeATs  Array of this particular molecule's atom types.
	 * @param moleculeNPs  The needed params of this molecule.
	 * @param isDepad      Whether to depad the output.
	 * @return Flattened 2d deriv array.
	 */
	public static double[] combine(double[][] derivs, int[] atomTypes,
								   int[][] neededParams, int[] moleculeATs,
								   int[][] moleculeNPs, boolean isDepad) {
		if (isDepad) derivs = depad(derivs, moleculeNPs);

		double[][] paddedDerivs = new double[atomTypes.length][];
		for (int i = 0; i < atomTypes.length; i++) {
			for (int j = 0; j < moleculeATs.length; j++) {
				if (atomTypes[i] == moleculeATs[j]) {
					paddedDerivs[i] = derivs[j];
				}
				else {
					paddedDerivs[i] = new double[neededParams[i].length];
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

	/**
	 * Fills up all gradient matrices, will be multithreaded if experimental
	 * geometry gradient computations are required as those take lots of time
	 * for they cannot be computed analytically.
	 * @return this
	 */
	public ParamGradient compute() {
		totalGradients =
				new double[s.getRm().mats.length][Solution.maxParamNum];
		if (analytical && (datum[1] != 0 || datum[2] != 0))
			computeBatchedDerivs(0, 0);

		if (!analytical || isExpAvail) {
			List<RecursiveAction> subtasks = new ArrayList<>();

			for (int Z = 0; Z < s.getRm().mats.length; Z++) {
				for (int paramNum : s.getRm().mnps[Z]) {
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
			for (int Z = 0; Z < s.getRm().mats.length; Z++) {
				for (int paramNum : s.getRm().mnps[Z]) {
					computeGradient(Z, paramNum);
				}
			}
		}

		return this;
	}

	/**
	 * Computes the gradients for one atom index or paramNum either
	 * analytically or by finite difference. Will also compute geometry
	 * derivatives if available.
	 *
	 * @param ZI       Atom index.
	 * @param paramNum Param number.
	 */
	protected void computeGradient(int ZI, int paramNum) {
		System.out.println("computing: " + ZI + " " + paramNum);
		Solution sPrime = null;
		if (!analytical) {
			sPrime = constructSPrime(ZI, paramNum);
		}

		computeHFDeriv(ZI, paramNum, sPrime);
		if (datum[1] != 0 && datum[2] != 0) {
			computeDipoleDeriv(ZI, paramNum, true, sPrime);
			computeIEDeriv(ZI, paramNum, sPrime);
		}
		else if (datum[1] != 0) {
			computeDipoleDeriv(ZI, paramNum, true, sPrime);
		}
		else if (datum[2] != 0) {
			computeDipoleDeriv(ZI, paramNum, false, sPrime);
			computeIEDeriv(ZI, paramNum, sPrime);
		}
		System.out.println("computing2: " + ZI + " " + paramNum);

		if (isExpAvail) computeGeomDeriv(ZI, paramNum);
	}

	/**
	 * Computes geometry derivative for one parameter of one atom.
	 *
	 * @param ZI       Atom index.
	 * @param paramNum Param number.
	 */
	protected void computeGeomDeriv(int ZI, int paramNum) {
		Solution sExpPrime = constructSExpPrime(ZI, paramNum);
		System.out.println("computing3: " + ZI + " " + paramNum);
		double sum = 0;
		for (int i = 0; i < sExpPrime.atoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				double d = findGrad(sExpPrime, i, j);
				sum += d * d;
			}
		}
		double geomGradient = 627.5 * Math.sqrt(sum);
		geomDerivs[ZI][paramNum] =
				1 / LAMBDA * (geomGradient - e.geomGradient);
		totalGradients[ZI][paramNum] +=
				0.000098 * e.geomGradient * geomDerivs[ZI][paramNum];
	}

	/**
	 * Computes error function of primary Solution object.
	 */
	protected void initializeErrorFunction() {
		if (datum[1] != 0) e.addDipoleError(datum[1]);
		if (datum[2] != 0) e.addIEError(datum[2]);

		if (this.sExp != null) {
			isExpAvail = true;
			geomDerivs = new double[s.getRm().mats.length]
					[Solution.maxParamNum];
			e.createExpGeom(this.sExp);
			e.addGeomError();
		}
	}

	/**
	 * Initializes result arrays as needed, so as to cut down on
	 * initialization time.
	 */
	protected void initializeArrays() {
		int atomLength = s.getRm().mats.length;
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

	/**
	 * Compiles all necessary fock matrices into one array before using the
	 * Pople algorithm, for faster computation. This function is the only
	 * thing that's not computed on a Z, paramNum level. Will not compute
	 * anything before the firstZIndex and the firstParamNum. Has not been
	 * implemented for unrestricted yet.
	 *
	 * @param firstZIndex   First atom index to compute, inclusive.
	 * @param firstParamNum First param number to compute, inclusive.
	 */
	protected abstract void computeBatchedDerivs(int firstZIndex,
												 int firstParamNum);

	/**
	 * Computes heat of formation derivative for one parameter of one atom.
	 *
	 * @param ZI       Atom index.
	 * @param paramNum Param number.
	 * @param sPrime   Perturbed s, used for finite difference if `analytical`
	 *                 is false.
	 */
	protected abstract void computeHFDeriv(int ZI, int paramNum,
										   Solution sPrime);

	/**
	 * Computes dipole derivative for one parameter of one atom. Also computes
	 * Hf derivs.
	 *
	 * @param ZI       Atom index.
	 * @param paramNum Param number.
	 * @param full     Sets whether dipole derivatives themselves are computed.
	 * @param sPrime   Perturbed s, used for finite difference if `analytical`
	 *                 is false.
	 */
	protected abstract void computeDipoleDeriv(int ZI, int paramNum,
											   boolean full,
											   Solution sPrime);

	/**
	 * Computes ionization energy derivative for one parameter of one atom.
	 * Requires computeDipoleDeriv to be called previously.
	 *
	 * @param ZI       Atom index.
	 * @param paramNum Param number.
	 * @param sPrime   Perturbed s, used for finite difference if `analytical`
	 *                 is false.
	 */
	protected abstract void computeIEDeriv(int ZI, int paramNum,
										   Solution sPrime);

	/**
	 * @param ZI       Atom index.
	 * @param paramNum Param number.
	 * @return Perturbed s for a certain param number of a certain atom index.
	 */
	protected abstract Solution constructSPrime(int ZI, int paramNum);

	/**
	 * @param ZI       Atom index.
	 * @param paramNum Param number.
	 * @return Perturbed sExp for a certain param number of a certain atom
	 * index.
	 */
	protected abstract Solution constructSExpPrime(int ZI, int paramNum);

	/**
	 * @param sExpPrime Perturbed sExp object.
	 * @param i         Index of atom, NOT unique atoms.
	 * @param xyz       Coordinate to find gradient of.
	 * @return Geometry gradient for that coordinate of that atom.
	 */
	protected abstract double findGrad(Solution sExpPrime, int i, int xyz);

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
