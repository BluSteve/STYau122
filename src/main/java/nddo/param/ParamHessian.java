package nddo.param;

import nddo.Constants;
import nddo.solution.Solution;
import org.apache.commons.lang3.ArrayUtils;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.stream.IntStream;

public class ParamHessian {
	protected final ParamGradient g;
	protected final Solution s, sExp;
	protected final double[] datum;
	protected final boolean restricted;
	protected double[][] hessian;
	protected boolean analytical;

	/**
	 * Constructs a ParamHessian object, not a time-intensive process.
	 *
	 * @param g The point of this class is to find the Hessian matrix
	 *          of the total error (that is, the errors of heat of
	 *          formation, dipole, ionization energy, and optionally
	 *          experimental geometry) with respect to various
	 *          parameters. g contains a Solution object which
	 *          represents what the program thinks is an optimal
	 *          solution to this molecule, and is not altered
	 *          throughout the process.
	 */
	private ParamHessian(ParamGradient g, boolean restricted) {
		this.g = g;
		s = g.s;
		datum = g.datum;
		sExp = g.sExp;
		analytical = g.analytical;
		this.restricted = restricted;
	}

	/**
	 * Takes in a pre-computed ParamGradient object.
	 *
	 * @param g ParamGradient object that contains the primary Solution
	 *          object.
	 * @return Uncomputed ParamHessian object, either restricted or not
	 * depending on type
	 * of g.
	 */
	public static ParamHessian from(ParamGradient g) {
		return new ParamHessian(g, g instanceof ParamGradientR);
	}

	/**
	 * Get Hessian matrix from basic ingredients.
	 *
	 * @param s     Primary Solution object.
	 * @param datum Reference data of Hf, dipole, I.E.
	 * @param sExp  Optional experimental geometry Solution object, may be
	 *              null.
	 * @return Uncomputed ParamHessian object, either restricted or not
	 * depending on type of s and sExp.
	 */
	public static ParamHessian of(Solution s, double[] datum, Solution sExp) {
		return ParamHessian.from(ParamGradient.of(s, datum, sExp).compute());
	}

	/**
	 * Gets flattened upper triangular matrix of
	 * Hessian.
	 *
	 * @param hessian 2dArray representing the
	 *                Hessian matrix. Should be
	 *                symmetrical.
	 * @return The upper triangular
	 * matrix flattened.
	 */
	public static double[] utify(double[][] hessian) {
		double[] hessianUT =
				new double[(hessian.length + 1) * hessian.length / 2];
		for (int i = 0; i < hessian.length; i++) {
			int p = 0;
			if (i >= 1) {
				p = i * (hessian.length * 2 - i + 1) / 2;
			}
			//noinspection ManualArrayCopy
			for (int j = i; j < hessian.length; j++)
				hessianUT[p + j - i] = hessian[i][j];
		}
		return hessianUT;
	}

	/**
	 * Returns hessian with all non-differentiated terms removed, padded for a
	 * particular training set.
	 *
	 * @param atomTypes    This is NOT the atom types of this particular
	 *                     molecule, but rather the atom types of the
	 *                     training set.
	 * @param neededParams This is NOT the neededParams of this particular
	 *                     molecule, but rather the neededParams of all the
	 *                     atoms in the training set.
	 * @return A Hessian matrix with all undifferentiated parameter elements
	 * (all of these are 0.0, but not all 0.0 elements are undifferentiated).
	 */
	public static double[][] padHessian(double[][] hessian, int[] mats,
										int[] atomTypes,
										int[][] neededParams) {
		List<List<Double>> padded = new ArrayList<>();
		int currentRowIndex = 0;
		for (int rowAT : atomTypes) {
			if (Utils.hasAtomType(mats, rowAT)) {
				for (int j = 0; j < Solution.maxParamNum; j++) {
					ArrayList<Double> row = new ArrayList<>();
					int currentColIndex = 0;
					for (int colAT : atomTypes) {
						if (Utils.hasAtomType(mats, colAT)) {
							for (int l = 0; l < Solution.maxParamNum; l++) {
								row.add(hessian[currentRowIndex][currentColIndex]);
								currentColIndex++;
							}
						}
						else {
							for (int l = 0; l < Solution.maxParamNum; l++) {
								row.add(0.0);
							}
						}
					}
					padded.add(row);
					currentRowIndex++;
				}
			}
			else {
				for (int j = 0; j < Solution.maxParamNum; j++) {
					padded.add(Arrays.asList(ArrayUtils.toObject(
							new double[Solution.maxParamNum *
									atomTypes.length])));
				}
			}
		}

		/*
		 * All parameters that are needed in the Hessian, stored linearly. I.e.
		 * can exceed MaxParamNum but not exceed MaxParamNum * number of
		 * different atoms.
		 */
		List<Integer> flattenedNPs = new ArrayList<>(Solution.maxParamNum *
				neededParams.length);
		int last = -1;
		for (int[] nps : neededParams) {
			for (int np : nps) {
				flattenedNPs.add(last + np + 1);
			}
			if (flattenedNPs.size() > 0)
				last = flattenedNPs.get(flattenedNPs.size() - 1);
		}
		double[][] unpadded =
				new double[flattenedNPs.size()][flattenedNPs.size()];
		int q = 0;
		for (int i : flattenedNPs) {
			int w = 0;
			for (int j : flattenedNPs) {
				unpadded[q][w] = padded.get(i).get(j);
				w++;
			}
			q++;
		}

		return unpadded;
	}

	/**
	 * Computes each row of the Hessian in a parallel manner. Very
	 * time-intensive.
	 *
	 * @return this
	 */
	public ParamHessian compute() {
		hessian = new double[s.rm.mnps.length * Solution.maxParamNum]
				[s.rm.mnps.length * Solution.maxParamNum];
		List<RecursiveAction> subtasks = new ArrayList<>();

		for (int ZIndex2 = 0; ZIndex2 < s.rm.mats.length; ZIndex2++) {
			for (int paramNum2 : s.rm.mnps[ZIndex2]) {
				int finalZIndex2 = ZIndex2;
				subtasks.add(new RecursiveAction() {
					@Override
					protected void compute() {
						computeRow(finalZIndex2, paramNum2);

					}
				});
			}
		}

		ForkJoinTask.invokeAll(subtasks);

		return this;
	}

	/**
	 * Computes one row of the Hessian in a parallel manner if each row would
	 * take a substantial amount of time to evaluate, i.e. analytical
	 * evaluation is turned off or experimental geometry is available.
	 *
	 * @param ZIndex2   The atom index of the row.
	 * @param paramNum2 The param number of the row, together with the atom
	 *                  index can determine the absolute row number.
	 */
	private void computeRow(int ZIndex2, int paramNum2) {
		ParamGradient gPrime = constructGPrime(ZIndex2, paramNum2);

		if (analytical && (datum[1] != 0 || datum[2] != 0))
			gPrime.computeBatchedDerivs(ZIndex2, paramNum2);

		List<RecursiveAction> subtasks = new ArrayList<>(s.rm.mats.length * Solution.maxParamNum);
		for (int ZIndex1 = ZIndex2; ZIndex1 < s.rm.mats.length; ZIndex1++) {
			for (int paramNum1 = paramNum2; paramNum1 < Solution.maxParamNum; paramNum1++) {
				boolean needed = false;
				for (int p : s.rm.mnps[ZIndex1]) {
					if (paramNum1 == p) {
						needed = true;
						break;
					}
				}

				if (needed) {
					if (!analytical || gPrime.isExpAvail) {
						int ZIndex = ZIndex1;
						int paramNum = paramNum1;
						subtasks.add(new RecursiveAction() {
							@Override
							protected void compute() {
								computeElement(gPrime, ZIndex2, paramNum2, ZIndex, paramNum);
							}
						});
					}
					else {
						computeElement(gPrime, ZIndex2, paramNum2, ZIndex1, paramNum1);
					}
				}
			}
		}
		if (subtasks.size() > 0) ForkJoinTask.invokeAll(subtasks);
	}

	private ParamGradient constructGPrime(int ZIndex, int paramNum) {
		return ParamGradient.of(s.withNewAtoms(Utils.perturbAtomParams(s.atoms,
						s.rm.mats[ZIndex], paramNum)), datum,
				sExp != null ? sExp.withNewAtoms(Utils.perturbAtomParams(sExp.atoms,
						sExp.rm.mats[ZIndex], paramNum)) : null);
	}

	/**
	 * Computes one element of the Hessian matrix.
	 *
	 * @param gPrime    The perturbed ParamGradient object used for this row
	 *                  of the Hessian
	 * @param ZIndex2   The row atom index.
	 * @param paramNum2 The row param number.
	 * @param ZIndex1   The column atom index.
	 * @param paramNum1 The column param number.
	 */
	private void computeElement(ParamGradient gPrime, int ZIndex2,
								int paramNum2, int ZIndex1, int paramNum1) {
		gPrime.computeGradient(ZIndex1, paramNum1);

		double element = (gPrime.getTotalGradients()[ZIndex1][paramNum1] -
				g.getTotalGradients()[ZIndex1][paramNum1]) / Constants.LAMBDA;

		int i2 = ZIndex2 * Solution.maxParamNum + paramNum2;
		int i1 = ZIndex1 * Solution.maxParamNum + paramNum1;

		hessian[i2][i1] = element;
		hessian[i1][i2] = element;
	}

	public double[][] getHessian() {
		return hessian;
	}

	public ParamGradient getG() {
		return g;
	}

	public ParamErrorFunction getE() {
		return g.getE();
	}

	@Deprecated
	public void computeSequentially() {
		hessian = new double[s.rm.mats.length * Solution.maxParamNum]
				[s.rm.mats.length * Solution.maxParamNum];

		// ZIndex2 and paramNum2 together denote the row number
		for (int ZIndex2 = 0; ZIndex2 < s.rm.mats.length; ZIndex2++) {
			for (int paramNum2 : s.rm.mnps[ZIndex2]) {
				computeRow(ZIndex2, paramNum2);
			}
		}
	}

	/**
	 * Gets the mse between analytical and finite difference.
	 */
	@Deprecated
	public double getAnalyticalError() {
		double sum = 0;
		try {
			this.compute();
			double[] a = ParamHessian.utify(this.hessian);
			double[][] b = new double[hessian.length][0];
			for (int i = 0; i < hessian.length; i++) b[i] = hessian[i].clone();

			analytical = !analytical;
			this.compute();
			sum = IntStream.range(0, a.length)
					.mapToDouble(i -> (a[i] -
							ParamHessian.utify(this.hessian)[i]) *
							(a[i] - ParamHessian.utify(this.hessian)[i]))
					.sum();
			analytical = !analytical;
			hessian = b;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return sum;
	}
}
