package nddo.param;

import nddo.math.PopleThiel;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Batcher;

import java.util.Arrays;

import static nddo.param.ParamSecondDerivative.*;

public class ParamHessianNew implements IParamHessian {
	final ParamGradientNew pg;
	final Solution s, sExp;
	final double[] datum;
	final boolean rhf, hasDip, hasIE, hasGeom;
	final int nAtomTypes, nParams;

	final double[][] hessian;
	final SimpleMatrix[][][][][] Fstatic2s, dD2statics, staticMatrices, PhiMatrices;

	public ParamHessianNew(ParamGradientNew pg) {
		this.pg = pg;

		s = pg.s;
		rhf = pg.rhf;
		SolutionR sr = rhf ? (SolutionR) s : null;
		SolutionU su = !rhf ? (SolutionU) s : null;

		sExp = pg.s;
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
		Fstatic2s = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
		dD2statics = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
		staticMatrices = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
		PhiMatrices = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];

		int[][] flat = new int[nanp * nanp][5]; // Z1, Z2, param1, param2;
		int[][] flatAll = new int[nanp * nanp][5]; // Z1, Z2, param1, param2;

		int i = 0;
		int iAll = 0;
		for (int ZI1 = 0; ZI1 < nAtomTypes; ZI1++) {
			for (int p1 : mnps[ZI1]) {
				for (int ZI2 = ZI1; ZI2 < nAtomTypes; ZI2++) {
					for (int p2 : mnps[ZI2]) {
						if (p2 >= p1) {
							flatAll[iAll][0] = ZI1;
							flatAll[iAll][1] = ZI2;
							flatAll[iAll][2] = p1;
							flatAll[iAll][3] = p2;
							flatAll[iAll][4] = iAll;
							iAll++;

							if (p1 != 0 && p2 != 0 && p1 != 7 && p2 != 7) {
								flat[i][0] = ZI1;
								flat[i][1] = ZI2;
								flat[i][2] = p1;
								flat[i][3] = p2;
								flat[i][4] = i;
								i++;
							}
						}
					}
				}
			}
		}

		System.out.println("i = " + i);
		System.out.println("iAll = " + iAll);

		flat = Arrays.copyOfRange(flat, 0, i);
		flatAll = Arrays.copyOfRange(flatAll, 0, iAll);

		SimpleMatrix[] ptInputsArr = new SimpleMatrix[flat.length];
		SimpleMatrix[] ptInputsArrBeta = !rhf ? new SimpleMatrix[flat.length] : null;

		Batcher.consume(flat, subset -> {
			for (int[] ints : subset) {
				int ZI1 = ints[0];
				int ZI2 = ints[1];
				int Z1 = mats[ZI1];
				int Z2 = mats[ZI2];
				int p1 = ints[2];
				int p2 = ints[3];

				SimpleMatrix Hderiv2 = Hderiv2(s, Z1, p1, Z2, p2);

				if (rhf) {
					SimpleMatrix[] Fstatic2 = Fstatic2s[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							Hderiv2.plusi(Gderiv2static(sr, Z1, p1, Z2, p2))
					};

					SimpleMatrix[] dD2static = dD2statics[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							densityderiv2static(sr, pg.xMatrices[ZI1][p1][0], pg.xMatrices[ZI2][p2][0])
					};

					SimpleMatrix[] PhiMatrix = PhiMatrices[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							staticFockDeriv(sr, Fstatic2[0],
									pg.densityDerivs[ZI1][p1][0], pg.densityDerivs[ZI2][p2][0],
									dD2static[0], Z1, p1, Z2, p2)
					};

					// dD2response precursor
					SimpleMatrix[] staticMatrix = staticMatrices[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							staticMatrix(sr, PhiMatrix[0], pg.FDerivs[ZI1][p1][0], pg.FDerivs[ZI2][p2][0],
									pg.xMatrices[ZI1][p1][0], pg.xMatrices[ZI2][p2][0])
					};

					ptInputsArr[ints[4]] = staticMatrix[0].extractMatrix(0, s.rm.nOccAlpha, s.rm.nOccAlpha,
							s.rm.nOrbitals);
				}
				else {
					SimpleMatrix[] Gderiv2static = Gderiv2static(su, Z1, p1, Z2, p2);

					Fstatic2s[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							Gderiv2static[0].plusi(Hderiv2), Gderiv2static[1].plusi(Hderiv2)
					};
				}
			}
		});

		SimpleMatrix[] dD2responses = rhf ? // flat.length
				Batcher.apply(ptInputsArr,
						subset -> {
							SimpleMatrix[] sms = PopleThiel.pople(sr, subset);
							SimpleMatrix[] results = new SimpleMatrix[sms.length];

							for (int j = 0; j < sms.length; j++) {
								results[j] = PopleThiel.densityDeriv(sr, sms[j]);
							}

							return results;
						}) : null;

		SimpleMatrix[][] dD2responsesU = !rhf ?
				Batcher.apply(new SimpleMatrix[][]{ptInputsArr, ptInputsArrBeta},
						subset -> {
							SimpleMatrix[] sms = PopleThiel.thiel(su, subset[0], subset[1]);
							SimpleMatrix[][] results = new SimpleMatrix[sms.length][];

							for (int j = 0; j < sms.length; j++) {
								results[j] = PopleThiel.densityDeriv(su, sms[j]);
							}

							return results;
						}) : null;

		int j = 0;
		for (int[] ints : flatAll) {
			int ZI1 = ints[0];
			int ZI2 = ints[1];
			int Z1 = mats[ZI1];
			int Z2 = mats[ZI2];
			int p1 = ints[2];
			int p2 = ints[3];

			if (Z1 == Z2 && p1 == 0 && p2 == 0) {
				double HfDeriv2 = alphaHfderiv2(s, Z1);

				addHfToHessian(ZI1, p1, ZI2, p2, HfDeriv2);
			}
			else if (p1 == 0 || p2 == 0 || p1 == 7 || p2 == 7) {
				addHfToHessian(ZI1, p1, ZI2, p2, 0);
			}
			else if (rhf) {
				SimpleMatrix dD2response = dD2responses[j++];
				SimpleMatrix densityDeriv2 = dD2response.plus(dD2statics[ZI1][ZI2][p1][p2][0]);

				double HfDeriv2 = MNDOHFDeriv2(sr, Z1, p1, Z2, p2,
						pg.staticDerivs[ZI1][0][p1], pg.staticDerivs[ZI1][1][p1], pg.densityDerivs[ZI2][p2][0], 0);

				addHfToHessian(ZI1, p1, ZI2, p2, HfDeriv2);


				if (hasDip) {
					double dipoleDeriv2 = MNDODipoleDeriv2(sr,
							pg.densityDerivs[ZI1][p1][0], pg.densityDerivs[ZI2][p2][0],
							densityDeriv2, Z1, p1, Z2, p2);

					addDipoleToHessian(ZI1, p1, ZI2, p2, dipoleDeriv2);
				}

				if (hasIE) {
					SimpleMatrix Phi = PhiMatrices[ZI1][ZI2][p1][p2][0];
					SimpleMatrix R = PopleThiel.responseMatrix(sr, dD2response);

					SimpleMatrix totalderiv = staticMatrices[ZI1][ZI2][p1][p2][0].plus(sr.Ct.mult(R).mult(sr.C));

					double IEDeriv2 = -homoDeriv2(sr,
							pg.xMatrices[ZI1][p1][0], pg.xMatrices[ZI2][p2][0], totalderiv,
							pg.FDerivs[ZI1][p1][0], pg.FDerivs[ZI2][p2][0],
							Phi.plus(R));

					addIEToHessian(ZI1, p1, ZI2, p2, IEDeriv2);
				}
			}
			else {

			}
		}
	}

	private void addToHessian(int ZI1, int p1, int ZI2, int p2, double x) {
		int i1 = ZI1 * nParams + p1;
		int i2 = ZI2 * nParams + p2;

		hessian[i1][i2] += x;
		if (i1 != i2) hessian[i2][i1] += x;
	}

	private void addHfToHessian(int ZI1, int p1, int ZI2, int p2, double x) {
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

	@Override
	public Solution getS() {
		return s;
	}

	@Override
	public double[][] getHessian() {
		return hessian;
	}
}
