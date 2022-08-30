package nddo.param;

import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import tools.Batcher;

import java.util.Arrays;

public class ParamHessianStewart implements IParamHessian {
	public final ParamGradientNew pg;
	final Solution s, sExp;
	final double[] datum;
	final boolean rhf, hasDip, hasIE, hasGeom;
	final int nAtomTypes, nParams;

	final double[][] hessian;
	private final SolutionR sr;
	private final SolutionU su;

	public ParamHessianStewart(Solution s, double[] datum, Solution sExp) {
		this(new ParamGradientNew(s, datum, sExp));
	}

	public ParamHessianStewart(ParamGradientNew pg) {
		this.pg = pg;

		s = pg.s;
		rhf = pg.rhf;
		sr = rhf ? (SolutionR) s : null;
		su = !rhf ? (SolutionU) s : null;

		sExp = pg.sExp;
		datum = pg.datum;
		hasDip = pg.hasDip;
		hasIE = pg.hasIE;
		hasGeom = pg.hasGeom;

		nAtomTypes = pg.nAtomTypes;
		nParams = pg.nParams;
		int nanp = nAtomTypes * nParams;
		int[] mats = s.rm.mats;
		int[][] mnps = s.rm.mnps;

		hessian = new double[nanp][nanp];
		int[][] flatAll = new int[nanp * nanp][6]; // Z1, Z2, param1, param2;

		int i = 0;
		int iAll = 0;
		for (int ZI1 = 0; ZI1 < nAtomTypes; ZI1++) {
			for (int p1 : mnps[ZI1]) {
				for (int ZI2 = ZI1; ZI2 < nAtomTypes; ZI2++) {
					for (int p2 : mnps[ZI2]) {
						flatAll[iAll][0] = ZI1;
						flatAll[iAll][1] = ZI2;
						flatAll[iAll][2] = p1;
						flatAll[iAll][3] = p2;
						flatAll[iAll][4] = i;
						flatAll[iAll][5] = iAll;
						iAll++;
					}
				}
			}
		}

		flatAll = Arrays.copyOfRange(flatAll, 0, iAll);

		Batcher.consume(flatAll, 1, subset -> {
			for (int[] ints : subset) {
				int ZI1 = ints[0];
				int ZI2 = ints[1];
				int Z1 = mats[ZI1];
				int Z2 = mats[ZI2];
				int p1 = ints[2];
				int p2 = ints[3];

				if (Z1 == Z2 && p1 == 0 && p2 == 0) {
					addHfToHessian(ZI1, p1, ZI2, p2, 0);
				}
				else if (p1 == 0 || p2 == 0 || p1 == 7 || p2 == 7) {
					addHfToHessian(ZI1, p1, ZI2, p2, 0);
				}
				else if (rhf) {
					addHfToHessian(ZI1, p1, ZI2, p2, 0);

					if (hasDip || hasIE) {
						if (hasDip) {
							addDipoleToHessian(ZI1, p1, ZI2, p2, 0);
						}

						if (hasIE) {
							addIEToHessian(ZI1, p1, ZI2, p2, 0);
						}
					}
				}
				else {
					addHfToHessian(ZI1, p1, ZI2, p2, 0);

					if (hasDip || hasIE) {
						if (hasDip) {
							addDipoleToHessian(ZI1, p1, ZI2, p2, 0);
						}

						if (hasIE) {
							addIEToHessian(ZI1, p1, ZI2, p2, 0);
						}
					}
				}
			}
		});

		if (hasGeom) {
			Batcher.consume(flatAll, 1, subset -> {
				for (int[] ints : subset) {
					s.rm.getLogger().trace("Starting batched derivs {} for ParamHessian geom", ints);

					int ZI1 = ints[0];
					int ZI2 = ints[1];
					int p1 = ints[2];
					int p2 = ints[3];

					if (p1 != 7 && p2 != 7) {
						addGeomToHessian(ZI1, p1, ZI2, p2, 0);
					}

					s.rm.getLogger().trace("Finished batched derivs {} for ParamHessian geom", ints);
				}
			});
		}
	}

	private void addToHessian(int ZI1, int p1, int ZI2, int p2, double x) {
		int i1 = ZI1 * nParams + p1;
		int i2 = ZI2 * nParams + p2;

		hessian[i1][i2] += x;
		if (i1 != i2) hessian[i2][i1] += x;
	}

	private void addHfToHessian(int ZI1, int p1, int ZI2, int p2, double x) { // x is 2nd derivative
		addToHessian(ZI1, p1, ZI2, p2, 2 * (pg.HfDerivs[ZI1][p1] * pg.HfDerivs[ZI2][p2] +
				(s.hf - datum[0]) * x));
	}

	private void addDipoleToHessian(int ZI1, int p1, int ZI2, int p2, double x) {
		addToHessian(ZI1, p1, ZI2, p2, 800 * (pg.dipoleDerivs[ZI1][p1] * pg.dipoleDerivs[ZI2][p2] +
				(s.dipole - datum[1]) * x));
	}

	private void addIEToHessian(int ZI1, int p1, int ZI2, int p2, double x) {
		addToHessian(ZI1, p1, ZI2, p2, 200 * (pg.IEDerivs[ZI1][p1] * pg.IEDerivs[ZI2][p2] -
				(s.homo + datum[2]) * x));
	}

	private void addGeomToHessian(int ZI1, int p1, int ZI2, int p2, double x) {
		addToHessian(ZI1, p1, ZI2, p2,
				0.25 * 2 * (pg.geomDerivs[ZI1][p1] * pg.geomDerivs[ZI2][p2] + pg.e.geomGradMag * x));
	}

	@Override
	public Solution getS() {
		return s;
	}

	@Override
	public ParamErrorFunction getE() {
		return pg.e;
	}

	@Override
	public double[][] getHessian() {
		return hessian;
	}
}
