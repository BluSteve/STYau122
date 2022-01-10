package nddo.param;

import nddo.Constants;
import nddo.State;
import nddo.math.PopleThiel;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Batcher;

import java.util.Arrays;

import static nddo.param.ParamGeometrySecondDerivative.gradDeriv2Alpha;
import static nddo.param.ParamGeometrySecondDerivative.gradderiv2;
import static nddo.param.ParamSecondDerivative.*;

public class ParamHessianNew implements IParamHessian {
	public final ParamGradientNew pg;
	final Solution s, sExp;
	final double[] datum;
	final boolean rhf, hasDip, hasIE, hasGeom;
	final int nAtomTypes, nParams;

	final double[][] hessian;
	final SimpleMatrix[][][][][] Fstatic2s, dD2statics, staticMatrices, PhiMatrices;
	private final SolutionR sr;
	private final SolutionU su;
	SimpleMatrix[][][][] gHVectorDerivs;
	SimpleMatrix[][][][][] Fstatic2sExp, dD2staticsExp, staticMatricesExp, PhiMatricesExp;

	public ParamHessianNew(Solution s, double[] datum, Solution sExp) {
		this(new ParamGradientNew(s, datum, sExp));
	}

	public ParamHessianNew(ParamGradientNew pg) {
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
		Fstatic2s = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
		dD2statics = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
		staticMatrices = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
		PhiMatrices = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];

		if (hasGeom) {
			Fstatic2sExp = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
			dD2staticsExp = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
			staticMatricesExp = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
			PhiMatricesExp = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];
			gHVectorDerivs = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams];
		}

		int[][] flat = new int[nanp * nanp][6]; // Z1, Z2, param1, param2;
		int[][] flatAll = new int[nanp * nanp][6]; // Z1, Z2, param1, param2;

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
							flatAll[iAll][4] = i;
							flatAll[iAll][5] = iAll;
							iAll++;

							if (p1 != 0 && p2 != 0 && p1 != 7 && p2 != 7) {
								flat[i][0] = ZI1;
								flat[i][1] = ZI2;
								flat[i][2] = p1;
								flat[i][3] = p2;
								flat[i][4] = i;
								flat[i][5] = iAll;
								i++;
							}
						}
					}
				}
			}
		}

		flat = Arrays.copyOfRange(flat, 0, i);
		flatAll = Arrays.copyOfRange(flatAll, 0, iAll);

		SimpleMatrix[] dD2responses = null;
		SimpleMatrix[][] dD2responsesU = null;

		if (hasIE || hasDip) {
			SimpleMatrix[] ptInputsArr = new SimpleMatrix[flat.length];
			SimpleMatrix[] ptInputsArrBeta = !rhf ? new SimpleMatrix[flat.length] : null;

			computeBatched(s, flat, ptInputsArr, ptInputsArrBeta, false);

			dD2responses = rhf ? computeDensityDerivs(sr, ptInputsArr) : null;

			dD2responsesU = !rhf ? computeDensityDerivs(su, ptInputsArr, ptInputsArrBeta) : null;
		}

		SimpleMatrix[] finalDD2responses = dD2responses;
		SimpleMatrix[][] finalDD2responsesU = dD2responsesU;
		Batcher.consume(flatAll, 1, subset -> {
			for (int[] ints : subset) {
				int ZI1 = ints[0];
				int ZI2 = ints[1];
				int Z1 = mats[ZI1];
				int Z2 = mats[ZI2];
				int p1 = ints[2];
				int p2 = ints[3];
				int j = ints[4];

				if (Z1 == Z2 && p1 == 0 && p2 == 0) {
					double HfDeriv2 = alphaHfderiv2(s, Z1);

					addHfToHessian(ZI1, p1, ZI2, p2, HfDeriv2);
				}
				else if (p1 == 0 || p2 == 0 || p1 == 7 || p2 == 7) {
					addHfToHessian(ZI1, p1, ZI2, p2, 0);
				}
				else if (rhf) {
					double HfDeriv2 = HfDeriv2(sr, Z1, p1, Z2, p2,
							pg.staticDerivs[ZI1][0][p1], pg.staticDerivs[ZI1][1][p1], pg.densityDerivs[ZI2][p2][0], 0);

					addHfToHessian(ZI1, p1, ZI2, p2, HfDeriv2);

					if (hasDip || hasIE) {
						SimpleMatrix dD2response = finalDD2responses[j];
						SimpleMatrix densityDeriv2 = dD2response.plus(dD2statics[ZI1][ZI2][p1][p2][0]);

						if (hasDip) {
							double dipoleDeriv2 = dipoleDeriv2(sr,
									pg.densityDerivs[ZI1][p1][0], pg.densityDerivs[ZI2][p2][0],
									densityDeriv2, Z1, p1, Z2, p2);

							addDipoleToHessian(ZI1, p1, ZI2, p2, dipoleDeriv2);
						}

						if (hasIE) {
							SimpleMatrix Phi = PhiMatrices[ZI1][ZI2][p1][p2][0];
							SimpleMatrix R = PopleThiel.responseMatrix(sr, dD2response);

							SimpleMatrix totalderiv =
									staticMatrices[ZI1][ZI2][p1][p2][0].plus(sr.Ct.mult(R).mult(sr.C));

							double IEDeriv2 = -homoDeriv2(sr,
									pg.xMatrices[ZI1][p1][0], pg.xMatrices[ZI2][p2][0], totalderiv,
									pg.FDerivs[ZI1][p1][0], pg.FDerivs[ZI2][p2][0],
									Phi.plus(R));

							addIEToHessian(ZI1, p1, ZI2, p2, IEDeriv2);
						}
					}
				}
				else {
					double HfDeriv2 = HfDeriv2(su, Z1, p1, Z2, p2,
							pg.staticDerivs[ZI1][0][p1], pg.staticDerivs[ZI1][1][p1], pg.staticDerivs[ZI1][2][p1],
							pg.densityDerivs[ZI2][p2][0], pg.densityDerivs[ZI2][p2][1], 0);

					addHfToHessian(ZI1, p1, ZI2, p2, HfDeriv2);

					if (hasDip || hasIE) {
						SimpleMatrix[] dD2response = finalDD2responsesU[j];
						SimpleMatrix[] densityDeriv2 = new SimpleMatrix[]{
								dD2response[0].plus(dD2statics[ZI1][ZI2][p1][p2][0]),
								dD2response[1].plus(dD2statics[ZI1][ZI2][p1][p2][1])
						};

						if (hasDip) {
							double dipoleDeriv2 = dipoleDeriv2(su,
									pg.densityDerivs[ZI1][p1][0].plus(pg.densityDerivs[ZI1][p1][1]),
									pg.densityDerivs[ZI2][p2][0].plus(pg.densityDerivs[ZI2][p2][1]),
									densityDeriv2[0].plus(densityDeriv2[1]), Z1, p1, Z2, p2);

							addDipoleToHessian(ZI1, p1, ZI2, p2, dipoleDeriv2);
						}

						if (hasIE) {
							SimpleMatrix[] Phi = PhiMatrices[ZI1][ZI2][p1][p2];
							SimpleMatrix[] R = PopleThiel.responseMatrix(su, dD2response);

							SimpleMatrix totalderiv =
									staticMatrices[ZI1][ZI2][p1][p2][0].plus(su.Cta.mult(R[0]).mult(su.Ca));

							double IEDeriv2 = -homoDeriv2(su,
									pg.xMatrices[ZI1][p1][0], pg.xMatrices[ZI2][p2][0], totalderiv,
									pg.FDerivs[ZI1][p1][0], pg.FDerivs[ZI2][p2][0],
									Phi[0].plus(R[0]));

							addIEToHessian(ZI1, p1, ZI2, p2, IEDeriv2);
						}
					}
				}
			}
		});

		if (hasGeom) {
			SimpleMatrix[] ptInputsArr = new SimpleMatrix[flat.length];
			SimpleMatrix[] ptInputsArrBeta = !rhf ? new SimpleMatrix[flat.length] : null;

			computeBatched(sExp, flat, ptInputsArr, ptInputsArrBeta, true);

			SimpleMatrix[] dD2responsesExp = rhf ? computeDensityDerivs((SolutionR) sExp, ptInputsArr) : null;

			SimpleMatrix[][] dD2responsesUExp =
					!rhf ? computeDensityDerivs((SolutionU) sExp, ptInputsArr, ptInputsArrBeta) : null;

			Batcher.consume(flatAll, 1, subset -> {
				for (int[] ints : subset) {
					s.rm.getLogger().trace("Starting batched derivs {} for ParamHessian geom", ints);

					int ZI1 = ints[0];
					int ZI2 = ints[1];
					int Z1 = mats[ZI1];
					int Z2 = mats[ZI2];
					int p1 = ints[2];
					int p2 = ints[3];
					int j = ints[4];

					SimpleMatrix deriv = gHVectorDerivs[ZI1][ZI2][p1][p2] = new SimpleMatrix(sExp.atoms.length * 3, 1);

					if (Z1 == Z2 && p1 == 0 && p2 == 0) {
						for (int atomNum = 0, k = 0; atomNum < sExp.atoms.length; atomNum++) {
							for (int tau = 0; tau < 3; tau++, k++) {
								deriv.set(k, gradDeriv2Alpha(sExp, atomNum, tau, mats[ZI1]));
							}
						}
					}
					else if (p1 != 0 && p1 != 7 && p2 != 0 && p2 != 7) {
						if (rhf) {
							SimpleMatrix densityDeriv2 = dD2responsesExp[j].plus(dD2staticsExp[ZI1][ZI2][p1][p2][0]);

							for (int atomNum = 0, k = 0; atomNum < sExp.atoms.length; atomNum++) {
								for (int tau = 0; tau < 3; tau++, k++) {
									deriv.set(k, gradderiv2((SolutionR) sExp, atomNum, tau, Z1, p1,
											Z2, p2, pg.densityDerivsExp[ZI1][p1][0], pg.densityDerivsExp[ZI2][p2][0],
											densityDeriv2));
								}
							}
						}
						else {
							SimpleMatrix[] densityDeriv2 = new SimpleMatrix[]{
									dD2responsesUExp[j][0].plus(dD2staticsExp[ZI1][ZI2][p1][p2][0]),
									dD2responsesUExp[j][1].plus(dD2staticsExp[ZI1][ZI2][p1][p2][1])
							};

							for (int atomNum = 0, k = 0; atomNum < sExp.atoms.length; atomNum++) {
								for (int tau = 0; tau < 3; tau++, k++) {
									deriv.set(k, gradderiv2((SolutionU) sExp, atomNum, tau, Z1, p1, Z2, p2,
											pg.densityDerivsExp[ZI1][p1][0], pg.densityDerivsExp[ZI1][p1][1],
											pg.densityDerivsExp[ZI2][p2][0], pg.densityDerivsExp[ZI2][p2][1],
											densityDeriv2[0], densityDeriv2[1]));
								}
							}
						}
					}

					if (p1 != 7 && p2 != 7) {
						addGeomToHessian(ZI1, p1, ZI2, p2, (Constants.KCAL *
								(deriv.dot(pg.e.geomGradVector) +
										pg.gGVectorDerivs[ZI1][p1].dot(pg.gGVectorDerivs[ZI2][p2])) -
								pg.geomDerivs[ZI1][p1] / Constants.KCAL * pg.geomDerivs[ZI2][p2]) /
								pg.e.geomGradMag * Constants.KCAL);
					}

					s.rm.getLogger().trace("Finished batched derivs {} for ParamHessian geom", ints);
				}
			});
		}
	}

	private static SimpleMatrix[] computeDensityDerivs(SolutionR sr, SimpleMatrix[] ptInputsArr) {
		return Batcher.apply(ptInputsArr, State.config.poplethiel_batch_size,
				subset -> {
					SimpleMatrix[] sms = PopleThiel.pt(sr, subset);
					SimpleMatrix[] results = new SimpleMatrix[sms.length];

					for (int j = 0; j < sms.length; j++) {
						results[j] = PopleThiel.densityDeriv(sr, sms[j]);
					}

					return results;
				});
	}

	private static SimpleMatrix[][] computeDensityDerivs(SolutionU su, SimpleMatrix[] ptInputsArr,
														 SimpleMatrix[] ptInputsArrBeta) {
		return Batcher.apply(ptInputsArr, ptInputsArrBeta, SimpleMatrix[][].class, State.config.poplethiel_batch_size,
				(subseta, subsetb) -> {
					SimpleMatrix[] sms = PopleThiel.pt(su, subseta, subsetb);
					SimpleMatrix[][] results = new SimpleMatrix[sms.length][];

					for (int j = 0; j < sms.length; j++) {
						results[j] = PopleThiel.densityDeriv(su, sms[j]);
					}

					return results;
				});
	}

	private void computeBatched(Solution s, int[][] flat, SimpleMatrix[] ptInputsArr, SimpleMatrix[] ptInputsArrBeta,
								boolean geom) {
		Batcher.consume(flat, 1, subset -> {
			for (int[] ints : subset) {
				final SolutionR sr = rhf ? (SolutionR) s : null;
				final SolutionU su = !rhf ? (SolutionU) s : null;

				int ZI1 = ints[0];
				int ZI2 = ints[1];
				int Z1 = s.rm.mats[ZI1];
				int Z2 = s.rm.mats[ZI2];
				int p1 = ints[2];
				int p2 = ints[3];

				SimpleMatrix[][][] densityDerivs, FDerivs, xMatrices;
				SimpleMatrix[][][][][] Fstatic2s, dD2statics, PhiMatrices, staticMatrices;
				if (geom) {
					densityDerivs = pg.densityDerivsExp;
					FDerivs = pg.FDerivsExp;
					xMatrices = pg.xMatricesExp;

					Fstatic2s = Fstatic2sExp;
					dD2statics = dD2staticsExp;
					PhiMatrices = PhiMatricesExp;
					staticMatrices = staticMatricesExp;
				}
				else {
					densityDerivs = pg.densityDerivs;
					FDerivs = pg.FDerivs;
					xMatrices = pg.xMatrices;

					Fstatic2s = this.Fstatic2s;
					dD2statics = this.dD2statics;
					PhiMatrices = this.PhiMatrices;
					staticMatrices = this.staticMatrices;
				}

				SimpleMatrix Hderiv2 = Hderiv2(s, Z1, p1, Z2, p2);

				if (rhf) {
					SimpleMatrix[] Fstatic2 = Fstatic2s[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							Gderiv2static(sr, Z1, p1, Z2, p2).plusi(Hderiv2)
					};

					SimpleMatrix[] dD2static = dD2statics[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							densityDeriv2static(sr, xMatrices[ZI1][p1][0], xMatrices[ZI2][p2][0])
					};

					SimpleMatrix[] PhiMatrix = PhiMatrices[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							staticFockDeriv(sr, Fstatic2[0],
									densityDerivs[ZI1][p1][0], densityDerivs[ZI2][p2][0],
									dD2static[0], Z1, p1, Z2, p2)
					};

					// dD2response precursor
					SimpleMatrix[] staticMatrix = staticMatrices[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							staticMatrix(sr, PhiMatrix[0], FDerivs[ZI1][p1][0], FDerivs[ZI2][p2][0],
									xMatrices[ZI1][p1][0], xMatrices[ZI2][p2][0])
					};

					ptInputsArr[ints[4]] =
							staticMatrix[0].extractMatrix(0, s.rm.nOccAlpha, s.rm.nOccAlpha, s.rm.nOrbitals);
				}
				else {
					SimpleMatrix[] Gderiv2static = Gderiv2static(su, Z1, p1, Z2, p2);

					SimpleMatrix[] Fstatic2 = Fstatic2s[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
							Gderiv2static[0].plusi(Hderiv2), Gderiv2static[1].plusi(Hderiv2)
					};

					SimpleMatrix[] dD2static = dD2statics[ZI1][ZI2][p1][p2] =
							densityDeriv2static(su, xMatrices[ZI1][p1], xMatrices[ZI2][p2]);

					SimpleMatrix[] PhiMatrix = PhiMatrices[ZI1][ZI2][p1][p2] =
							staticFockDeriv(su, Fstatic2,
									densityDerivs[ZI1][p1], densityDerivs[ZI2][p2],
									dD2static, Z1, p1, Z2, p2);

					SimpleMatrix[] staticMatrix = staticMatrices[ZI1][ZI2][p1][p2] =
							staticMatrix(su, PhiMatrix, FDerivs[ZI1][p1], FDerivs[ZI2][p2],
									xMatrices[ZI1][p1], xMatrices[ZI2][p2]);

					ptInputsArr[ints[4]] =
							staticMatrix[0].extractMatrix(0, s.rm.nOccAlpha, s.rm.nOccAlpha, s.rm.nOrbitals);

					ptInputsArrBeta[ints[4]] =
							staticMatrix[1].extractMatrix(0, s.rm.nOccBeta, s.rm.nOccBeta, s.rm.nOrbitals);
				}
			}
		});
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
