package nddo.param;

import nddo.Constants;
import nddo.State;
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
import static nddo.param.ParamGeometryDerivative.gradDeriv;
import static nddo.param.ParamGeometryDerivative.gradDerivAlpha;

public final class ParamGradientNew implements IParamGradient {
	final Solution s, sExp;
	final ParamErrorFunction e;
	final double[] datum;
	final boolean rhf, hasDip, hasIE, hasGeom;
	final int nAtomTypes, nParams;

	double[][] HfDerivs, dipoleDerivs, IEDerivs, geomDerivs, totalGradients;

	SimpleMatrix[][] xVectors;
	SimpleMatrix[][][] staticDerivs, densityDerivs, FDerivs, xMatrices;

	SimpleMatrix[][] xVectorsExp, gGVectorDerivs;
	SimpleMatrix[][][] staticDerivsExp, densityDerivsExp, FDerivsExp, xMatricesExp;

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

		staticDerivs = new SimpleMatrix[nAtomTypes][][];
		xVectors = new SimpleMatrix[nAtomTypes][];
		densityDerivs = new SimpleMatrix[nAtomTypes][nParams][];

		FDerivs = new SimpleMatrix[nAtomTypes][nParams][];
		xMatrices = new SimpleMatrix[nAtomTypes][nParams][];

		if (hasDip) {
			dipoleDerivs = new double[nAtomTypes][nParams];
			e.addDipoleError(datum[1]);
		}

		if (hasIE) {
			IEDerivs = new double[nAtomTypes][nParams];
			e.addIEError(datum[2]);
		}

		if (hasGeom) {
			e.createExpGeom(sExp);
			e.addGeomError();

			geomDerivs = new double[nAtomTypes][nParams];

			staticDerivsExp = new SimpleMatrix[nAtomTypes][][];
			xVectorsExp = new SimpleMatrix[nAtomTypes][];
			densityDerivsExp = new SimpleMatrix[nAtomTypes][nParams][];

			FDerivsExp = new SimpleMatrix[nAtomTypes][nParams][];
			xMatricesExp = new SimpleMatrix[nAtomTypes][nParams][];

			gGVectorDerivs = new SimpleMatrix[nAtomTypes][nParams];
		}


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

					if (rhf) {
						SimpleMatrix responseMatrix = PopleThiel.responseMatrix(sr, densityDerivs[ZI][paramNum][0]);

						SimpleMatrix Fao = responseMatrix.plusi(staticDerivs[ZI][1][paramNum]);
						FDerivs[ZI][paramNum] = new SimpleMatrix[]{Fao};

						SimpleMatrix Fmo = sr.Ct.mult(Fao).mult(sr.C);

						SimpleMatrix xMatrix = xMatrix(sr, Fmo);
						xMatrices[ZI][paramNum] = new SimpleMatrix[]{xMatrix};

						if (hasIE) IEDerivs[ZI][paramNum] = -homoDeriv(sr, xMatrix, Fmo);
					}
					else {
						SimpleMatrix[] responseMatrices = PopleThiel.responseMatrix(su, densityDerivs[ZI][paramNum]);

						SimpleMatrix Faoa = responseMatrices[0].plusi(staticDerivs[ZI][1][paramNum]);
						SimpleMatrix Faob = responseMatrices[1].plusi(staticDerivs[ZI][2][paramNum]);
						FDerivs[ZI][paramNum] = new SimpleMatrix[]{Faoa, Faob};

						SimpleMatrix Fmoa = su.Cta.mult(Faoa).mult(su.Ca);
						SimpleMatrix Fmob = su.Ctb.mult(Faob).mult(su.Cb);

						xMatrices[ZI][paramNum] = ParamDerivative.xMatrix(su, Fmoa, Fmob);

						if (hasIE) IEDerivs[ZI][paramNum] = -homoDeriv(su, xMatrices[ZI][paramNum][0], Fmoa);
					}

					if (hasIE) addIEGrad(ZI, paramNum);
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

			Batcher.consume(params.toArray(new int[0][0]), 1, subset -> {
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

							if (rhf) {
								SolutionR sExpr = (SolutionR) sExp;

								SimpleMatrix responseMatrix =
										PopleThiel.responseMatrix(sExpr, densityDerivsExp[ZI][paramNum][0]);

								SimpleMatrix Fao = responseMatrix.plusi(staticDerivsExp[ZI][1][paramNum]);
								FDerivsExp[ZI][paramNum] = new SimpleMatrix[]{Fao};

								SimpleMatrix Fmo = sExpr.Ct.mult(Fao).mult(sExpr.C);

								SimpleMatrix xMatrix = xMatrix(sExpr, Fmo);
								xMatricesExp[ZI][paramNum] = new SimpleMatrix[]{xMatrix};
							}
							else {
								SolutionU sExpu = (SolutionU) sExp;

								SimpleMatrix[] responseMatrices =
										PopleThiel.responseMatrix(sExpu, densityDerivsExp[ZI][paramNum]);

								SimpleMatrix Faoa = responseMatrices[0].plusi(staticDerivsExp[ZI][1][paramNum]);
								SimpleMatrix Faob = responseMatrices[1].plusi(staticDerivsExp[ZI][2][paramNum]);
								FDerivsExp[ZI][paramNum] = new SimpleMatrix[]{Faoa, Faob};

								SimpleMatrix Fmoa = sExpu.Cta.mult(Faoa).mult(sExpu.Ca);
								SimpleMatrix Fmob = sExpu.Ctb.mult(Faob).mult(sExpu.Cb);

								xMatricesExp[ZI][paramNum] = ParamDerivative.xMatrix(sExpu, Fmoa, Fmob);
							}
						}

						geomDerivs[ZI][paramNum] =
								Constants.KCAL * Constants.KCAL * deriv.dot(e.geomGradVector) / e.geomGradMag;
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

		for (int ZI = 0, i = 0; ZI < atomLength; ZI++) { // todo make granular
			sdTarget[ZI] = rhf ? staticDeriv(sr, s.rm.mats[ZI], 0) : staticDeriv(su, s.rm.mats[ZI], 0);

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
					Batcher.apply(aggFaUnpad, State.config.poplethiel_batch_size,
							subset -> PopleThiel.pt(sr, PopleThiel.toMO(sr.CtOcc, sr.CVirt, subset))) :
					Batcher.apply(aggFaUnpad, aggFbUnpad, SimpleMatrix[].class, State.config.poplethiel_batch_size,
							(a, b) -> PopleThiel.pt(su, PopleThiel.toMO(su.CtaOcc, su.CaVirt, a),
									PopleThiel.toMO(su.CtbOcc, su.CbVirt, b)));

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
		totalGradients[ZI][paramNum] += 0.25 * 2 * e.geomGradMag * geomDerivs[ZI][paramNum];
	}

	public ParamErrorFunction getE() {
		return e;
	}

	@Override
	public Solution getS() {
		return s;
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
