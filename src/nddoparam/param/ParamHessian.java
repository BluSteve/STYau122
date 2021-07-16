package nddoparam.param;

import nddoparam.Solution;
import scf.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;

public abstract class ParamHessian {
	protected ParamGradient g;
	protected Solution s, sExp;
	protected int[] atomTypes;
	protected double[] datum;
	protected double[][] hessian;
	protected boolean analytical;

	public ParamHessian(Solution s, double[] datum,
						Solution sExp, boolean analytical, int[] atomTypes) {
		this.s = s;
		this.datum = datum;
		this.sExp = sExp;
		this.analytical = analytical;
		this.atomTypes = atomTypes;
	}

	public static double[] getHessianUT(double[][] unpadded) {
		double[] hessianUT =
				new double[(unpadded.length + 1) * unpadded.length / 2];
		for (int i = 0; i < unpadded.length; i++) {
			int p = 0;
			if (i >= 1) {
				p = i * (unpadded.length * 2 - i + 1) / 2;
			}
			for (int j = i; j < unpadded.length; j++)
				hessianUT[p + j - i] = unpadded[i][j];
		}
		return hessianUT;
	}

	public void computeHessian() {
		hessian = new double[atomTypes.length * Solution.maxParamNum]
				[atomTypes.length * Solution.maxParamNum];
		List<int[]> ZandPNs = new ArrayList<>(hessian.length);

		for (int ZIndex2 = 0; ZIndex2 < s.getUniqueZs().length; ZIndex2++) {
			for (int paramNum2 : s.getNeededParams()
					[s.getUniqueZs()[ZIndex2]]) {
				ZandPNs.add(new int[]{ZIndex2, paramNum2});
			}
		}

		int cores = Runtime.getRuntime().availableProcessors();
		ForkJoinPool pool =
				new ForkJoinPool(Utils.findTightestTriplet(cores)[1]);

		pool.submit(() -> ZandPNs.parallelStream()
				.forEach(request -> {
					int ZIndex2 = request[0];
					int paramNum2 = request[1];
					computeHessianRow(ZIndex2, paramNum2);
				}));
	}

	protected void computeHessianRow(int ZIndex2, int paramNum2) {
		ParamGradient gPrime = constructGPrime(ZIndex2, paramNum2);

		if (analytical && (datum[1] != 0 || datum[2] != 0))
			gPrime.computeBatchedDerivs(ZIndex2, paramNum2);

		boolean needed;
		for (int ZIndex1 = ZIndex2; ZIndex1 < s.getUniqueZs().length;
			 ZIndex1++) {
			for (int paramNum1 = paramNum2;
				 paramNum1 < Solution.maxParamNum;
				 paramNum1++) {
				needed = false;
				for (int p : s.getNeededParams()[s
						.getUniqueZs()[ZIndex1]]) {
					if (paramNum1 == p) {
						needed = true;
						break;
					}
				}

				if (needed) {
					gPrime.computeGradient(ZIndex1, paramNum1);
					hessian[ZIndex2 * Solution.maxParamNum + paramNum2]
							[ZIndex1 * Solution.maxParamNum +
							paramNum1] =
							(gPrime.getTotalGradients()
									[ZIndex1][paramNum1] -
									g.getTotalGradients()
											[ZIndex1][paramNum1]) /
									Utils.LAMBDA;

					hessian[ZIndex1 * Solution.maxParamNum + paramNum1]
							[ZIndex2 * Solution.maxParamNum +
							paramNum2] =
							hessian[ZIndex2 * Solution.maxParamNum +
									paramNum2]
									[ZIndex1 * Solution.maxParamNum +
									paramNum1];
				}
			}
		}
	}

	public void computeHessianSequential() {
		hessian = new double[atomTypes.length * Solution.maxParamNum]
				[atomTypes.length * Solution.maxParamNum];

		for (int ZIndex2 = 0; ZIndex2 < s.getUniqueZs().length; ZIndex2++) {
			for (int paramNum2 : s.getNeededParams()[s
					.getUniqueZs()[ZIndex2]]) {
				computeHessianRow(ZIndex2, paramNum2);
			}
		}
	}

	protected abstract ParamGradient constructGPrime(int ZIndex, int paramNum);

	// returns hessian with all non-differentiated terms removed.
	// length is Main's paramLength ** 2.
	public double[][] getHessianUnpadded(int[] atomTypes,
										 int[][] neededParams,
										 int paramLength) {
		ArrayList<Integer> pList = new ArrayList<>();
		int last = -1;
		for (int[] i : neededParams) {
			for (int i1 : i) {
				pList.add(last + i1 + 1);
			}
			if (pList.size() > 0) last = pList.get(pList.size() - 1);
		}
		double[][] unpadded = new double[pList.size()][pList.size()];
		int q = 0;
		for (int i : pList) {
			int w = 0;
			for (int j : pList) {
				unpadded[q][w] = hessian[i][j];
				w++;
			}
			q++;
		}

		return unpadded;
	}

	public double[] getHessianUT() {
		return getHessianUT(getHessianUnpadded(null, null, 0));
	}

	public double[][] getHessian() {
		return hessian;
	}

	// Gets the mse between analytical and finite difference.
	public double getAnalyticalError() {
		double sum = 0;
		try {
			this.computeHessian();
			double[] a = getHessianUT().clone();
			double[][] b = new double[hessian.length][0];
			for (int i = 0; i < hessian.length; i++) b[i] = hessian[i].clone();

			analytical = !analytical;
			this.computeHessian();
			sum = IntStream.range(0, a.length)
					.mapToDouble(i -> (a[i] - getHessianUT()[i]) *
							(a[i] - getHessianUT()[i]))
					.sum();
			analytical = !analytical;
			hessian = b;
		} catch (Exception e) {
		}
		return sum;
	}

	public ParamErrorFunction getE() {
		return g.getE();
	}

	public ParamGradient getG() {
		return g;
	}

	public boolean isAnalytical() {
		return analytical;
	}

	public void setAnalytical(boolean analytical) {
		this.analytical = analytical;
	}
}
