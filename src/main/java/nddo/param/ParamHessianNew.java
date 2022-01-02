package nddo.param;

import nddo.math.PopleThiel;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Batcher;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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

		int[][] flat = new int[nanp * nanp][4]; // Z1, Z2, param1, param2;

		int i = 0;
		for (int ZI1 = 0; ZI1 < nAtomTypes; ZI1++) {
			for (int p1 : mnps[ZI1]) {
				if (p1 != 0 && p1 != 7) {
					for (int ZI2 = ZI1; ZI2 < nAtomTypes; ZI2++) {
						for (int p2 : mnps[ZI2]) {
							if (p2 != 0 && p2 != 7 && p2 >= p1) {
								flat[i][0] = ZI1;
								flat[i][1] = ZI2;
								flat[i][2] = p1;
								flat[i][3] = p2;
								i++;
							}
						}
					}
				}
			}
		}

		System.out.println("i = " + i);

		flat = Arrays.copyOfRange(flat, 0, i);

		List<SimpleMatrix> ptInputs = new ArrayList<>();
		List<SimpleMatrix> ptInputsBeta = !rhf ? new ArrayList<>() : null;

		for (int[] ints : flat) {
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

//				System.out.println(Arrays.toString(ints) + Fstatic2[0].elementSum());

				SimpleMatrix[] dD2static = dD2statics[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
						densityderiv2static(sr, pg.xMatrices[ZI1][p1][0], pg.xMatrices[ZI2][p2][0])
				};

				SimpleMatrix[] PhiMatrix = PhiMatrices[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
						staticFockDeriv(sr, Fstatic2[0],
								pg.densityDerivs[ZI1][p1][0], pg.densityDerivs[ZI2][p2][0],
								dD2static[0], Z1, p1, Z2, p2)
				};

//				SimpleMatrix[] staticMatrix = staticMatrices[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
//						staticMatrix(sr, PhiMatrix[0], pg.FDerivs[ZI1][p1][0], pg.FDerivs[ZI2][p2][0],
//								pg.xMatrices[ZI1][p1][0], pg.xMatrices[ZI2][p2][0])
//				};

				SimpleMatrix[] staticMatrix = staticMatrices[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
						staticMatrix(sr, Fstatic2[0], pg.staticDerivs[ZI1][0][p1],pg.staticDerivs[ZI2][0][p2], pg.xVectors[ZI1][p1], pg.xVectors[ZI2][p2], Z1, p1, Z2, p2)
				};


				ptInputs.add(staticMatrix[0].extractMatrix(0, s.rm.nOccAlpha, s.rm.nOccAlpha, s.rm.nOrbitals));

				//System.out.println(Arrays.toString(ints) + Fstatic2[0].elementSum());

			}
			else {
				SimpleMatrix[] Gderiv2static = Gderiv2static(su, Z1, p1, Z2, p2);

				Fstatic2s[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
						Gderiv2static[0].plusi(Hderiv2), Gderiv2static[1].plusi(Hderiv2)
				};
			}
		}

		SimpleMatrix[] dD2responses = rhf ? // flat.length
				Batcher.apply(ptInputs.toArray(new SimpleMatrix[0]),
						subset -> {
							SimpleMatrix[] sms = PopleThiel.pople(sr, subset);
							SimpleMatrix[] results = new SimpleMatrix[sms.length];

							for (int j = 0; j < sms.length; j++) {
								results[j] = PopleThiel.densityDeriv(sr, sms[j]);
							}

							return results;
						}) : null;

		SimpleMatrix[][] dD2responsesU = !rhf ?
				Batcher.apply(new SimpleMatrix[][]{ptInputs.toArray(new SimpleMatrix[0]),
								ptInputsBeta.toArray(new SimpleMatrix[0])},
						subset -> {
							SimpleMatrix[] sms = PopleThiel.thiel(su, subset[0], subset[1]);
							SimpleMatrix[][] results = new SimpleMatrix[sms.length][];

							for (int j = 0; j < sms.length; j++) {
								results[j] = PopleThiel.densityDeriv(su, sms[j]);
							}

							return results;
						}) : null;

		for (int j = 0; j < flat.length; j++) {
			int[] ints = flat[j];
			int ZI1 = ints[0];
			int ZI2 = ints[1];
			int Z1 = mats[ZI1];
			int Z2 = mats[ZI2];
			int p1 = ints[2];
			int p2 = ints[3];

			if (rhf) {
				SimpleMatrix dD2response = dD2responses[j];
				SimpleMatrix densityDeriv2 = dD2response.plusi(dD2statics[ZI1][ZI2][p1][p2][0]);
				SimpleMatrix Phi = PhiMatrices[ZI1][ZI2][p1][p2][0];

				double HfDeriv2 = MNDOHFDeriv2(sr, Z1, p1, Z2, p2,
						pg.staticDerivs[ZI1][0][p1], pg.staticDerivs[ZI1][1][p1], pg.densityDerivs[ZI2][p2][0], 0);

				addToHessian(ZI1, p1, ZI2, p2,
						2 * (pg.HfDerivs[ZI1][p1] * pg.HfDerivs[ZI2][p2] + (s.hf - datum[0]) * HfDeriv2));

				if (hasDip) {
					double dipoleDeriv2 = MNDODipoleDeriv2(sr,
							pg.densityDerivs[ZI1][p1][0], pg.densityDerivs[ZI2][p2][0],
									densityDeriv2, Z1, p1, Z2, p2);

					addToHessian(ZI1, p1, ZI2, p2, 800 * (pg.dipoleDerivs[ZI1][p1] * pg.dipoleDerivs[ZI2][p2] +
									(s.dipole - datum[1]) * dipoleDeriv2));
				}

				SimpleMatrix R = PopleThiel.responseMatrix(sr, dD2response);

				SimpleMatrix totalderivs = sr.Ct.mult(Phi.plus(R)).mult(sr.C);

//	e IEderiv2 = MNDO			do

			}
			else {

			}
		}

//		System.out.println("Fstatic = " + Arrays.deepToString(Fstatic));
	}

	private void addToHessian(int ZI1, int p1, int ZI2, int p2, double x) {
		int i1 = ZI1 * nParams + p1;
		int i2 = ZI2 * nParams + p2;

		hessian[i1][i2] += x;
		if (i1 != i2) hessian[i2][i1] += x;
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
