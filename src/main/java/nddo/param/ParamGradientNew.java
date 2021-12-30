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

public class ParamGradientNew implements IParamGradient {
	protected final Solution s, sExp;
	protected final ParamErrorFunction e;
	protected final double[] datum;
	protected final boolean rhf, hasDip, hasIE, hasGeom;
	protected final int atomLength, paramLength;

	protected double[][] HFDerivs, dipoleDerivs, IEDerivs, geomDerivs, totalGradients;
	protected SimpleMatrix[][] xVectors;
	protected SimpleMatrix[][][] densityDerivs, staticDerivs, FDerivs, xMatrices;

	public ParamGradientNew(Solution s, double[] datum, Solution sExp) {
		this.s = s;
		this.datum = datum;
		this.sExp = sExp;
		this.rhf = s instanceof SolutionR;

		atomLength = s.rm.mats.length;
		paramLength = Solution.maxParamNum;

		hasDip = datum[1] != 0;
		hasIE = datum[2] != 0;
		hasGeom = sExp != null;


		e = ParamErrorFunction.of(s, datum[0]);

		totalGradients = new double[atomLength][paramLength];
		HFDerivs = new double[atomLength][paramLength];

		if (hasDip || hasIE) {
			staticDerivs = new SimpleMatrix[atomLength][2][paramLength];
			xVectors = new SimpleMatrix[atomLength][paramLength];
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
			geomDerivs = new double[s.rm.mats.length][Solution.maxParamNum];
			e.createExpGeom(sExp);
			e.addGeomError();
		}


		if (hasDip || hasIE) {
			SolutionR sr = rhf ? (SolutionR) s : null;
			SolutionU su = !rhf ? (SolutionU) s : null;

			SimpleMatrix[] aggFa = new SimpleMatrix[atomLength * paramLength];
			SimpleMatrix[] aggFb = new SimpleMatrix[atomLength * paramLength];

			for (int ZI = 0, i = 0; ZI < atomLength; ZI++) {
				staticDerivs[ZI] = rhf ? MNDOStaticMatrixDeriv(sr, s.rm.mats[ZI], 0)
						: MNDOStaticMatrixDeriv(su, s.rm.mats[ZI], 0);

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
				for (int paramNum : s.rm.mnps[ZI]) {
					if (paramNum == 0 || paramNum == 7) {
						HFDerivs[ZI][paramNum] = HFDeriv(s, s.rm.mats[ZI], paramNum);
						addHFDeriv(ZI, paramNum);
					}
					else if (staticDerivs[ZI][0][paramNum] != null) {
						HFDerivs[ZI][paramNum] = rhf ?
								MNDOHFDeriv(sr, staticDerivs[ZI][0][paramNum], staticDerivs[ZI][1][paramNum]) :
								MNDOHFDeriv(su, staticDerivs[ZI][0][paramNum], staticDerivs[ZI][1][paramNum],
										staticDerivs[ZI][2][paramNum]);
						addHFDeriv(ZI, paramNum);

						if (hasDip || hasIE) {
							SimpleMatrix dd;
							if (rhf) {
								densityDerivs[ZI][paramNum] =
										new SimpleMatrix[]{PopleThiel.densityDeriv(sr, xVectors[ZI][paramNum])};
								dd = densityDerivs[ZI][paramNum][0];
							}
							else {
								densityDerivs[ZI][paramNum] = PopleThiel.densityDeriv(su, xVectors[ZI][paramNum]);
								dd = densityDerivs[ZI][paramNum][0].plus(densityDerivs[ZI][paramNum][1]);
							}

							if (hasDip) {
								dipoleDerivs[ZI][paramNum] = MNDODipoleDeriv(s, dd, s.rm.mats[ZI], paramNum);
								addDipoleDeriv(ZI, paramNum);
							}
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

							addIEDeriv(ZI, paramNum);
						}
					}
				}
			}
		}
		else {
			for (int ZI = 0; ZI < atomLength; ZI++) {
				for (int paramNum : s.rm.mnps[ZI]) {
					HFDerivs[ZI][paramNum] = HFDeriv(s, s.rm.mats[ZI], paramNum);
					addHFDeriv(ZI, paramNum);
				}
			}
		}
	}

	private void addHFDeriv(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
	}

	private void addDipoleDeriv(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 800 * (s.dipole - datum[1]) * dipoleDerivs[ZI][paramNum];
	}

	private void addIEDeriv(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 200 * -(s.homo + datum[2]) * IEDerivs[ZI][paramNum];
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
