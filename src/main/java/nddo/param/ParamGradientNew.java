package nddo.param;

import nddo.math.PopleThiel;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Batcher;
import tools.Utils;

import java.util.Arrays;

import static nddo.param.ParamDerivative.*;
import static nddo.param.ParamGeometryDerivative.gradderiv;

public class ParamGradientNew implements IParamGradient {
	protected final Solution s, sExp;
	protected final ParamErrorFunction e;
	protected final double[] datum;
	protected final boolean rhf, hasDip, hasIE, hasGeom;
	protected final int atomLength, paramLength;

	protected double[][] HFDerivs, dipoleDerivs, IEDerivs, geomDerivs, totalGradients;
	protected SimpleMatrix[][] xVectors, gGVectorDerivs;
	protected SimpleMatrix[][][] staticDerivs, densityDerivs, FDerivs, xMatrices;

	public ParamGradientNew(Solution s, double[] datum, Solution sExp) {
		this.s = s;
		this.datum = datum;
		this.sExp = sExp;
		this.rhf = s instanceof SolutionR;

		int[] mats = s.rm.mats;
		int[][] mnps = s.rm.mnps;
		atomLength = mats.length;
		paramLength = Solution.maxParamNum;

		hasDip = datum[1] != 0;
		hasIE = datum[2] != 0;
		hasGeom = sExp != null;


		e = ParamErrorFunction.of(s, datum[0]);

		totalGradients = new double[atomLength][paramLength];
		HFDerivs = new double[atomLength][paramLength];

		if (hasDip || hasIE || hasGeom) {
			staticDerivs = new SimpleMatrix[atomLength][][];
			xVectors = new SimpleMatrix[atomLength][];
			densityDerivs = new SimpleMatrix[atomLength][paramLength][];

			if (hasDip) {
				dipoleDerivs = new double[atomLength][paramLength];
				e.addDipoleError(datum[1]);
			}

			if (hasIE) {
				FDerivs = new SimpleMatrix[atomLength][paramLength][];
				xMatrices = new SimpleMatrix[atomLength][paramLength][];
				IEDerivs = new double[atomLength][paramLength];
				e.addIEError(datum[2]);
			}
		}

		if (hasGeom) {
			geomDerivs = new double[atomLength][paramLength];
			gGVectorDerivs = new SimpleMatrix[atomLength][paramLength];

			e.createExpGeom(sExp);
			e.addGeomError();
		}


		if (hasDip || hasIE || hasGeom) {
			SolutionR sr = rhf ? (SolutionR) s : null;
			SolutionU su = !rhf ? (SolutionU) s : null;

			SimpleMatrix[] aggFa = new SimpleMatrix[atomLength * paramLength];
			SimpleMatrix[] aggFb = new SimpleMatrix[atomLength * paramLength];

			for (int ZI = 0, i = 0; ZI < atomLength; ZI++) {
				staticDerivs[ZI] = rhf ? MNDOStaticMatrixDeriv(sr, mats[ZI], 0)
						: MNDOStaticMatrixDeriv(su, mats[ZI], 0);

				for (int j = 0; j < staticDerivs[ZI][1].length; j++) {
					aggFa[i] = staticDerivs[ZI][1][j];
					if (!rhf) aggFb[i] = staticDerivs[ZI][2][j];
					i++;
				}
			}

			SimpleMatrix[] aggFaUnpad = new SimpleMatrix[Utils.numNotNull(aggFa)];
			SimpleMatrix[] aggFbUnpad = !rhf ? new SimpleMatrix[Utils.numNotNull(aggFb)] : null;
			for (int i = 0, ui = 0; i < aggFa.length; i++) {
				if (aggFa[i] != null) {
					aggFaUnpad[ui] = aggFa[i];
					if (!rhf) aggFbUnpad[ui] = aggFb[i];
					ui++;
				}
			}

			if (aggFa.length > 0) {
				SimpleMatrix[] aggXUnpad = rhf ?
						Batcher.apply(aggFaUnpad,
								subset -> PopleThiel.pople(sr, PopleThiel.aoToMo(sr.CtOcc, sr.CVirt, subset))) :
						Batcher.apply(new SimpleMatrix[][]{aggFaUnpad, aggFbUnpad},
								subset -> PopleThiel.thiel(su, PopleThiel.aoToMo(su.CtaOcc, su.CaVirt, subset[0]),
										PopleThiel.aoToMo(su.CtbOcc, su.CbVirt, subset[1])));

				SimpleMatrix[] aggX = new SimpleMatrix[aggFa.length];

				for (int i = 0, ui = 0; i < aggFa.length; i++) {
					if (aggFa[i] != null) {
						aggX[i] = aggXUnpad[ui];
						ui++;
					}
				}

				for (int i = 0; i < atomLength; i++) {
					xVectors[i] = Arrays.copyOfRange(aggX, i * paramLength, (i + 1) * paramLength);
				}
			}


			for (int ZI = 0; ZI < atomLength; ZI++) {
				for (int paramNum : mnps[ZI]) {
					if (paramNum == 0 || paramNum == 7) {
						HFDerivs[ZI][paramNum] = HFDeriv(s, mats[ZI], paramNum);
						addHFGrad(ZI, paramNum);
					}
					else if (staticDerivs[ZI][0][paramNum] != null) {
						HFDerivs[ZI][paramNum] = rhf ?
								MNDOHFDeriv(sr, staticDerivs[ZI][0][paramNum], staticDerivs[ZI][1][paramNum]) :
								MNDOHFDeriv(su, staticDerivs[ZI][0][paramNum], staticDerivs[ZI][1][paramNum],
										staticDerivs[ZI][2][paramNum]);
						addHFGrad(ZI, paramNum);


						SimpleMatrix combinedDd;
						if (rhf) {
							densityDerivs[ZI][paramNum] =
									new SimpleMatrix[]{PopleThiel.densityDeriv(sr, xVectors[ZI][paramNum])};

							combinedDd = densityDerivs[ZI][paramNum][0];
						}
						else {
							densityDerivs[ZI][paramNum] = PopleThiel.densityDeriv(su, xVectors[ZI][paramNum]);

							combinedDd = densityDerivs[ZI][paramNum][0].plus(densityDerivs[ZI][paramNum][1]);
						}


						if (hasDip) {
							dipoleDerivs[ZI][paramNum] = MNDODipoleDeriv(s, combinedDd, mats[ZI], paramNum);
							addDipoleGrad(ZI, paramNum);
						}


						if (hasIE) {
							if (rhf) {
								SimpleMatrix responseMatrix =
										PopleThiel.responseMatrix(sr, densityDerivs[ZI][paramNum][0]);

								SimpleMatrix plus = staticDerivs[ZI][1][paramNum].plus(responseMatrix);
								SimpleMatrix F = sr.Ct.mult(plus).mult(sr.C);
								FDerivs[ZI][paramNum] = new SimpleMatrix[]{F};

								SimpleMatrix xMatrix = xmatrix(F, sr);
								xMatrices[ZI][paramNum] = new SimpleMatrix[]{xMatrix};

								IEDerivs[ZI][paramNum] = -MNDOHomoDerivNew(sr, xMatrix, plus);
							}
							else {
								SimpleMatrix[] responseMatrices =
										PopleThiel.responseMatrices(su, densityDerivs[ZI][paramNum]);

								SimpleMatrix plus = staticDerivs[ZI][1][paramNum].plus(responseMatrices[0]);
								SimpleMatrix Fa = su.Cta.mult(plus).mult(su.Ca);
								SimpleMatrix Fb = su.Ctb.mult(staticDerivs[ZI][2][paramNum].plus(responseMatrices[1]))
										.mult(su.Cb);
								FDerivs[ZI][paramNum] = new SimpleMatrix[]{Fa, Fb};

								xMatrices[ZI][paramNum] = xmatrices(Fa, Fb, su);

								IEDerivs[ZI][paramNum] = -MNDOHomoDerivNew(su, xMatrices[ZI][paramNum][0], plus);
							}

							addIEGrad(ZI, paramNum);
						}


						if (hasGeom) {
							SimpleMatrix deriv = gGVectorDerivs[ZI][paramNum] = new SimpleMatrix(atomLength * 3, 1);
							for (int atomnum = 0, i = 0; atomnum < atomLength; atomnum++) {
								for (int tau = 0; tau < 3; tau++) {
									deriv.set(i++, rhf ?
											gradderiv((SolutionR) sExp, atomnum, tau, mats[ZI], paramNum,
													densityDerivs[ZI][paramNum][0]) :
											gradderiv((SolutionU) sExp, atomnum, tau, mats[ZI], paramNum,
													densityDerivs[ZI][paramNum][0], densityDerivs[ZI][paramNum][1]));
								}
							}

							geomDerivs[ZI][paramNum] = 627.5 * 627.5 * deriv.dot(e.geomGradVector) / e.geomGradient;
							addGeomGrad(ZI, paramNum);
						}
					}
				}
			}
		}
		else {
			for (int ZI = 0; ZI < atomLength; ZI++) {
				for (int paramNum : mnps[ZI]) {
					HFDerivs[ZI][paramNum] = HFDeriv(s, mats[ZI], paramNum);
					addHFGrad(ZI, paramNum);
				}
			}
		}
	}

	private void addHFGrad(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
	}

	private void addDipoleGrad(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 800 * (s.dipole - datum[1]) * dipoleDerivs[ZI][paramNum];
	}

	private void addIEGrad(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 200 * -(s.homo + datum[2]) * IEDerivs[ZI][paramNum];
	}

	private void addGeomGrad(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 0.000098 * e.geomGradient * geomDerivs[ZI][paramNum];
	}

	@Override
	public Solution getS() {
		return s;
	}

	@Override
	public ParamErrorFunction getE() {
		return e;
	}

	@Override
	public double[][] getHFDerivs() { // todo uncapitalize F
		return HFDerivs;
	}

	@Override
	public double[][] getDipoleDerivs() {
		return dipoleDerivs;
	}

	@Override
	public double[][] getIEDerivs() {
		return IEDerivs;
	}

	@Override
	public double[][] getGeomDerivs() {
		return geomDerivs;
	}

	@Override
	public double[][] getTotalGradients() {
		return totalGradients;
	}
}
