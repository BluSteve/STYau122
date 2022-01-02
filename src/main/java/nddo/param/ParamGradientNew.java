package nddo.param;

import nddo.math.PopleThiel;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Batcher;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static nddo.param.ParamDerivative.*;
import static nddo.param.ParamGeometryDerivative.*;

public class ParamGradientNew implements IParamGradient {
	protected final Solution s, sExp;
	protected final ParamErrorFunction e;
	protected final double[] datum;
	protected final boolean rhf, hasDip, hasIE, hasGeom;
	protected final int nAtomTypes, nParams;

	protected double[][] HfDerivs, dipoleDerivs, IEDerivs, geomDerivs, totalGradients;

	protected SimpleMatrix[][] xVectors;
	protected SimpleMatrix[][][] staticDerivs, densityDerivs, FDerivs, xMatrices;

	protected SimpleMatrix[][] xVectorsExp, gGVectorDerivs;
	protected SimpleMatrix[][][] staticDerivsExp, densityDerivsExp;

	public ParamGradientNew(Solution s, double[] datum, Solution sExp) {
		this.s = s;
		this.datum = datum;
		this.sExp = sExp;
		this.rhf = s instanceof SolutionR;

		int[] mats = s.rm.mats;
		int[][] mnps = s.rm.mnps;
		nAtomTypes = mats.length; // how many unique atom types there are
		nParams = Solution.maxParamNum;

		hasDip = datum[1] != 0;
		hasIE = datum[2] != 0;
		hasGeom = sExp != null;

		e = ParamErrorFunction.of(s, datum[0]);

		totalGradients = new double[nAtomTypes][nParams];
		HfDerivs = new double[nAtomTypes][nParams];

		if (hasDip || hasIE || hasGeom) {
			staticDerivs = new SimpleMatrix[nAtomTypes][][];
			xVectors = new SimpleMatrix[nAtomTypes][];
			densityDerivs = new SimpleMatrix[nAtomTypes][nParams][];

			if (hasDip) {
				dipoleDerivs = new double[nAtomTypes][nParams];
				e.addDipoleError(datum[1]);
			}

			if (hasIE) {
				FDerivs = new SimpleMatrix[nAtomTypes][nParams][];
				xMatrices = new SimpleMatrix[nAtomTypes][nParams][];
				IEDerivs = new double[nAtomTypes][nParams];
				e.addIEError(datum[2]);
			}
		}

		if (hasGeom) {
			e.createExpGeom(sExp);
			e.addGeomError();

			geomDerivs = new double[nAtomTypes][nParams];

			staticDerivsExp = new SimpleMatrix[nAtomTypes][][];
			xVectorsExp = new SimpleMatrix[nAtomTypes][];
			densityDerivsExp = new SimpleMatrix[nAtomTypes][nParams][];

			gGVectorDerivs = new SimpleMatrix[nAtomTypes][nParams];
		}


		if (hasDip || hasIE) {
			SolutionR sr = rhf ? (SolutionR) s : null;
			SolutionU su = !rhf ? (SolutionU) s : null;

			computeBatchedDerivs(s, staticDerivs, xVectors);

			for (int ZI = 0; ZI < nAtomTypes; ZI++) {
				for (int paramNum : mnps[ZI]) {
					if (paramNum == 0 || paramNum == 7) {
						HfDerivs[ZI][paramNum] = ParamDerivative.HfDeriv(s, mats[ZI], paramNum);
						addHfGrad(ZI, paramNum);
					}
					else if (staticDerivs[ZI][0][paramNum] != null) {
						HfDerivs[ZI][paramNum] = rhf ?
								HfDeriv(sr, staticDerivs[ZI][0][paramNum], staticDerivs[ZI][1][paramNum]) :
								HfDeriv(su, staticDerivs[ZI][0][paramNum], staticDerivs[ZI][1][paramNum],
										staticDerivs[ZI][2][paramNum]);
						addHfGrad(ZI, paramNum);

						computeDensityDerivs(s, xVectors, ZI, paramNum, densityDerivs);

						if (hasDip) {
							SimpleMatrix combinedDd = rhf ? densityDerivs[ZI][paramNum][0] :
									densityDerivs[ZI][paramNum][0].plus(densityDerivs[ZI][paramNum][1]);

							dipoleDerivs[ZI][paramNum] = nddoDipoleDeriv(s, combinedDd, mats[ZI], paramNum);
							addDipoleGrad(ZI, paramNum);
						}

						if (hasIE) {
							if (rhf) {
								SimpleMatrix responseMatrix =
										PopleThiel.responseMatrix(sr, densityDerivs[ZI][paramNum][0]);

								SimpleMatrix plus = staticDerivs[ZI][1][paramNum].plus(responseMatrix);
								SimpleMatrix F = sr.Ct.mult(plus).mult(sr.C);
								FDerivs[ZI][paramNum] = new SimpleMatrix[]{F};

								SimpleMatrix xMatrix = xMatrix(sr, F);
								xMatrices[ZI][paramNum] = new SimpleMatrix[]{xMatrix};

								IEDerivs[ZI][paramNum] = -homoDeriv(sr, xMatrix, plus);
							}
							else {
								SimpleMatrix[] responseMatrices =
										PopleThiel.responseMatrices(su, densityDerivs[ZI][paramNum]);

								SimpleMatrix plus = staticDerivs[ZI][1][paramNum].plus(responseMatrices[0]);
								SimpleMatrix Fa = su.Cta.mult(plus).mult(su.Ca);
								SimpleMatrix Fb = su.Ctb.mult(staticDerivs[ZI][2][paramNum].plus(responseMatrices[1]))
										.mult(su.Cb);
								FDerivs[ZI][paramNum] = new SimpleMatrix[]{Fa, Fb};

								xMatrices[ZI][paramNum] = ParamDerivative.xMatrix(su, Fa, Fb);

								IEDerivs[ZI][paramNum] = -homoDeriv(su, xMatrices[ZI][paramNum][0], plus);
							}

							addIEGrad(ZI, paramNum);
						}
					}
				}
			}
		}
		else {
			for (int ZI = 0; ZI < nAtomTypes; ZI++) {
				for (int paramNum : mnps[ZI]) {
					HfDerivs[ZI][paramNum] = ParamDerivative.HfDeriv(s, mats[ZI], paramNum);
					addHfGrad(ZI, paramNum);
				}
			}
		}

		if (hasGeom) {
			computeBatchedDerivs(sExp, staticDerivsExp, xVectorsExp);

			List<int[]> params = new ArrayList<>(nAtomTypes * nParams);
			for (int ZI = 0; ZI < nAtomTypes; ZI++) {
				for (int j = 0; j < mnps[ZI].length; j++) {
					int paramNum = mnps[ZI][j];
					int[] ints = new int[]{ZI, paramNum};
					params.add(ints);
				}
			}

			Batcher.consume(params.toArray(new int[0][0]), subset -> {
				for (int[] ints : subset) {
					int ZI = ints[0];
					int paramNum = ints[1];
					if (paramNum == 0 || staticDerivsExp[ZI][0][paramNum] != null) {
						SimpleMatrix deriv = gGVectorDerivs[ZI][paramNum] = new SimpleMatrix(sExp.atoms.length * 3, 1);

						if (paramNum == 0) {
							for (int atomNum = 0, i = 0; atomNum < sExp.atoms.length; atomNum++) {
								for (int tau = 0; tau < 3; tau++, i++) {
									deriv.set(i, gradDerivAlpha(sExp, atomNum, tau, mats[ZI]));
								}
							}
						}
						else {
							computeDensityDerivs(sExp, xVectorsExp, ZI, paramNum, densityDerivsExp);

							for (int atomNum = 0, i = 0; atomNum < sExp.atoms.length; atomNum++) {
								for (int tau = 0; tau < 3; tau++, i++) {
									deriv.set(i, rhf ?
											gradDeriv((SolutionR) sExp, atomNum, tau, mats[ZI], paramNum,
													densityDerivsExp[ZI][paramNum][0]) :
											gradDeriv((SolutionU) sExp, atomNum, tau, mats[ZI], paramNum,
													densityDerivsExp[ZI][paramNum][0],
													densityDerivsExp[ZI][paramNum][1])
									);
								}
							}
						}

						geomDerivs[ZI][paramNum] = 627.5 * 627.5 * deriv.dot(e.geomGradVector) / e.geomGradMag;
						addGeomGrad(ZI, paramNum);
					}
				}
			});
		}
	}

	private static void computeDensityDerivs(Solution s, SimpleMatrix[][] xVectors, int ZI, int paramNum,
											 SimpleMatrix[][][] dDTarget) {
		dDTarget[ZI][paramNum] = s instanceof SolutionR ?
				new SimpleMatrix[]{PopleThiel.densityDeriv((SolutionR) s, xVectors[ZI][paramNum])} :
				PopleThiel.densityDeriv((SolutionU) s, xVectors[ZI][paramNum]);
	}

	private static void computeBatchedDerivs(Solution s, SimpleMatrix[][][] sdTarget, SimpleMatrix[][] xVTarget) {
		final int atomLength = s.rm.mats.length;
		final int paramLength = Solution.maxParamNum;
		final boolean rhf = s instanceof SolutionR;
		final SolutionR sr = rhf ? (SolutionR) s : null;
		final SolutionU su = !rhf ? (SolutionU) s : null;

		SimpleMatrix[] aggFa = new SimpleMatrix[atomLength * paramLength];
		SimpleMatrix[] aggFb = new SimpleMatrix[atomLength * paramLength];

		for (int ZI = 0, i = 0; ZI < atomLength; ZI++) {
			sdTarget[ZI] = rhf ? staticDeriv(sr, s.rm.mats[ZI], 0)
					: staticDeriv(su, s.rm.mats[ZI], 0);

			for (int j = 0; j < sdTarget[ZI][1].length; j++) {
				aggFa[i] = sdTarget[ZI][1][j];
				if (!rhf) aggFb[i] = sdTarget[ZI][2][j];
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
							subset -> PopleThiel.pople(sr, PopleThiel.toMO(sr.CtOcc, sr.CVirt, subset))) :
					Batcher.apply(new SimpleMatrix[][]{aggFaUnpad, aggFbUnpad},
							subset -> PopleThiel.thiel(su, PopleThiel.toMO(su.CtaOcc, su.CaVirt, subset[0]),
									PopleThiel.toMO(su.CtbOcc, su.CbVirt, subset[1])));

			SimpleMatrix[] aggX = new SimpleMatrix[aggFa.length];

			for (int i = 0, ui = 0; i < aggFa.length; i++) {
				if (aggFa[i] != null) {
					aggX[i] = aggXUnpad[ui];
					ui++;
				}
			}

			for (int i = 0; i < atomLength; i++) {
				xVTarget[i] = Arrays.copyOfRange(aggX, i * paramLength, (i + 1) * paramLength);
			}
		}
	}

	private void addHfGrad(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 2 * (s.hf - datum[0]) * HfDerivs[ZI][paramNum];
	}

	private void addDipoleGrad(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 800 * (s.dipole - datum[1]) * dipoleDerivs[ZI][paramNum];
	}

	private void addIEGrad(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 200 * -(s.homo + datum[2]) * IEDerivs[ZI][paramNum];
	}

	private void addGeomGrad(int ZI, int paramNum) {
		totalGradients[ZI][paramNum] += 0.000098 * e.geomGradMag * geomDerivs[ZI][paramNum];
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
	public double[][] getHfDerivs() {
		return HfDerivs;
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
