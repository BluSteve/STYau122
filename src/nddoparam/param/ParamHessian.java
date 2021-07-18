package nddoparam.param;

import nddoparam.Solution;
import org.apache.commons.lang3.ArrayUtils;
import scf.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.stream.IntStream;

public abstract class ParamHessian {
	protected ParamGradient g;
	protected Solution s, sExp;
	protected double[] datum;
	protected double[][] hessian;
	protected boolean analytical;

	/**
	 * Constructs a ParamHessian object, not a time-intensive process.
	 *
	 * @param s          The point of this class is to find the Hessian matrix
	 *                   of the total error (that is, the errors of heat of
	 *                   formation, dipole, ionization energy, and optionally
	 *                   experimental geometry) with respect to various
	 *                   parameters. This Solution object represents what the
	 *                   program thinks is an optimal solution to this
	 *                   molecule, and is not altered throughout the process.
	 * @param datum      Reference data using which we can derive the errors
	 *                   of heat of formation, dipole, and I.E.
	 * @param sExp       Solution created using the experimental geometry of
	 *                   the molecule. Used to compute the geometry error
	 *                   and gradients.
	 * @param analytical Boolean indicating whether analytical derivatives
	 *                   should be used when available. Should be true by
	 *                   default.
	 */
	public ParamHessian(Solution s, double[] datum, Solution sExp,
						boolean analytical) {
		this.s = s;
		this.datum = datum;
		this.sExp = sExp;
		this.analytical = analytical;
	}

	/**
	 * Gets flattened upper triangular matrix of Hessian.
	 *
	 * @param hessian 2dArray representing the Hessian matrix. Should be
	 *                symmetrical.
	 * @return The upper triangular matrix flattened.
	 */
	public static double[] getHessianUT(double[][] hessian) {
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
	 * Computes each row of the Hessian in a parallel manner.
	 */
	public void compute() {
		hessian = new double[s.getUniqueZs().length * Solution.maxParamNum]
				[s.getUniqueZs().length * Solution.maxParamNum];
		List<RecursiveAction> subtasks = new ArrayList<>();

		for (int ZIndex2 = 0; ZIndex2 < s.getUniqueZs().length; ZIndex2++) {
			for (int paramNum2 : s.getNeededParams()
					[s.getUniqueZs()[ZIndex2]]) {
				int finalZIndex2 = ZIndex2;
				subtasks.add(new RecursiveAction() {
					@Override
					protected void compute() {
						computeHessianRow(finalZIndex2, paramNum2);
					}
				});
			}
		}

		ForkJoinTask.invokeAll(subtasks);
	}

	/**
	 * Computes one row of the Hessian in a parallel manner the row has enough
	 * elements. Note that Hessian matrices are symmetrical, so if too few
	 * elements need to be computed, the overhead from multithreading becomes
	 * significant and sequential computation is preferred.
	 *
	 * @param ZIndex2   The atom index of the row.
	 * @param paramNum2 The param number of the row, together with the atom
	 *                  index can determine the absolute row number.
	 */
	protected void computeHessianRow(int ZIndex2, int paramNum2) {
		ParamGradient gPrime = constructGPrime(ZIndex2, paramNum2);

		if (analytical && (datum[1] != 0 || datum[2] != 0))
			gPrime.computeBatchedDerivs(ZIndex2, paramNum2);

		List<RecursiveAction> subtasks = new ArrayList<>();
		for (int ZIndex1 = ZIndex2; ZIndex1 < s.getUniqueZs().length;
			 ZIndex1++) {
			List<int[]> ZandPNs = new ArrayList<>();
			for (int paramNum1 = paramNum2; paramNum1 < Solution.maxParamNum;
				 paramNum1++) {
				boolean needed = false;
				for (int p : s.getNeededParams()[s.getUniqueZs()[ZIndex1]]) {
					if (paramNum1 == p) {
						needed = true;
						break;
					}
				}

				if (needed) ZandPNs.add(new int[]{ZIndex1, paramNum1});
			}

			// Multithread only if there is expGeom and at least 5 gradient
			// computations are needed
			if (sExp != null && ZandPNs.size() >= 5) {
				for (int[] ZandPN : ZandPNs) {
					int ZIndex = ZandPN[0];
					int paramNum = ZandPN[1];
					subtasks.add(new RecursiveAction() {
						@Override
						protected void compute() {
							gPrime.computeGradient(ZIndex,
									paramNum);
							int i2 =
									ZIndex2 * Solution.maxParamNum + paramNum2;
							int i1 = ZIndex * Solution.maxParamNum +
									paramNum;

							hessian[i2][i1] =
									(gPrime.getTotalGradients()
											[ZIndex]
											[paramNum] -
											g.getTotalGradients()
													[ZIndex]
													[paramNum])
											/ Utils.LAMBDA;
							hessian[i1][i2] = hessian[i2][i1];
						}
					});
				}
			}
			else {
				for (int[] ZandPN : ZandPNs) {
					gPrime.computeGradient(ZandPN[0], ZandPN[1]);
					int i2 = ZIndex2 * Solution.maxParamNum + paramNum2;
					int i1 = ZandPN[0] * Solution.maxParamNum + ZandPN[1];

					hessian[i2][i1] =
							(gPrime.getTotalGradients()[ZandPN[0]][ZandPN[1]] -
									g.getTotalGradients()[ZandPN[0]][ZandPN[1]])
									/ Utils.LAMBDA;
					hessian[i1][i2] = hessian[i2][i1];
				}
			}
		}
		if (subtasks.size() > 1) ForkJoinTask.invokeAll(subtasks);
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
	public double[][] getHessian(int[] atomTypes, int[][] neededParams) {
		List<List<Double>> padded = new ArrayList<>();
		int currentRowIndex = 0;
		for (int rowAT : atomTypes) {
			if (Utils.hasAtomType(s, rowAT)) {
				for (int j = 0; j < Solution.maxParamNum; j++) {
					ArrayList<Double> row = new ArrayList<>();
					int currentColIndex = 0;
					for (int colAT : atomTypes) {
						if (Utils.hasAtomType(s, colAT)) {
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

	public double[][] getHessianRaw() {
		return hessian;
	}

	public ParamGradient getG() {
		return g;
	}

	public ParamErrorFunction getE() {
		return g.getE();
	}

	protected abstract ParamGradient constructGPrime(int ZIndex, int paramNum);

	@Deprecated
	public void computeSequentially() {
		hessian = new double[s.getUniqueZs().length * Solution.maxParamNum]
				[s.getUniqueZs().length * Solution.maxParamNum];

		// ZIndex2 and paramNum2 together denote the row number
		for (int ZIndex2 = 0; ZIndex2 < s.getUniqueZs().length; ZIndex2++) {
			for (int paramNum2 : s.getNeededParams()[s
					.getUniqueZs()[ZIndex2]]) {
				computeHessianRow(ZIndex2, paramNum2);
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
			double[] a = ParamHessian.getHessianUT(this.hessian);
			double[][] b = new double[hessian.length][0];
			for (int i = 0; i < hessian.length; i++) b[i] = hessian[i].clone();

			analytical = !analytical;
			this.compute();
			sum = IntStream.range(0, a.length)
					.mapToDouble(i -> (a[i] -
							ParamHessian.getHessianUT(this.hessian)[i]) *
							(a[i] - ParamHessian.getHessianUT(this.hessian)[i]))
					.sum();
			analytical = !analytical;
			hessian = b;
		} catch (Exception ignored) {
		}
		return sum;
	}
}
