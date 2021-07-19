package nddoparam.param;

import nddoparam.GeometryDerivative;
import nddoparam.Solution;
import nddoparam.SolutionR;
import org.jblas.DoubleMatrix;
import scf.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

class ParamGradientR extends ParamGradient {
	protected ParamGradientR(SolutionR s, double[] datum, SolutionR sExp,
						  boolean analytical) {
		super(s, datum, sExp, analytical);
	}

	@Override
	protected void computeBatchedDerivs(int firstZIndex, int firstParamNum) {
		ArrayList<DoubleMatrix> aggregate =
				new ArrayList<>(s.getRm().mats.length * Solution.maxParamNum);
		for (int ZI = 0; ZI < s.getRm().mats.length; ZI++) {
			if (ZI == firstZIndex)
				staticDerivs[ZI] = ParamDerivative
						.MNDOStaticMatrixDeriv((SolutionR) s,
								s.getRm().mats[ZI],
								firstParamNum);
			else if (ZI < firstZIndex) {
				staticDerivs[ZI] =
						new DoubleMatrix[][]{
								new DoubleMatrix[Solution.maxParamNum],
								new DoubleMatrix[Solution.maxParamNum]};
			}
			else staticDerivs[ZI] = ParamDerivative
						.MNDOStaticMatrixDeriv((SolutionR) s,
								s.getRm().mats[ZI], 0);
			Collections.addAll(aggregate, staticDerivs[ZI][1]);
		}

		DoubleMatrix[] aggregateArray = new DoubleMatrix[aggregate.size()];
		for (int x = 0; x < aggregate.size(); x++) {
			aggregateArray[x] = aggregate.get(x);
		}
		DoubleMatrix[] xLimitedAggregate =
				ParamDerivative.xArrayLimitedPople((SolutionR) s,
						aggregateArray);
		int i = 0;
		for (int Z = 0; Z < s.getRm().mats.length; Z++) {
			xLimited[Z] =
					Arrays.copyOfRange(xLimitedAggregate,
							i * Solution.maxParamNum,
							i * Solution.maxParamNum + Solution.maxParamNum);
			i++;
		}
	}

	@Override
	protected void computeHFDeriv(int ZI, int paramNum, Solution sPrime) {
		if (analytical)
			HFDerivs[ZI][paramNum] = ParamDerivative
					.HFDeriv((SolutionR) s, s.getRm().mats[ZI], paramNum);
		else {
			assert sPrime != null;
			HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
		}
		totalGradients[ZI][paramNum] +=
				2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
	}

	@Override
	protected void computeDipoleDeriv(int ZI, int paramNum, boolean full,
									  Solution sPrime) {
		if (analytical) {
			// alpha and eisol have to be computed apart from the rest of
			// HFDerivs. I.e not as a part of dipole.
			if (paramNum == 0 || paramNum == 7) {
				HFDerivs[ZI][paramNum] = ParamDerivative
						.HFDeriv((SolutionR) s, s.getRm().mats[ZI],
								paramNum);
			}
			else if (staticDerivs[ZI][0][paramNum] != null ||
					staticDerivs[ZI][1][paramNum] != null) {
				HFDerivs[ZI][paramNum] =
						ParamDerivative.MNDOHFDeriv((SolutionR) s,
								staticDerivs[ZI][0][paramNum],
								staticDerivs[ZI][1][paramNum]);
				densityDerivs[ZI][paramNum] = ParamDerivative
						.densityDerivativeLimited((SolutionR) s,
								xLimited[ZI][paramNum]);
				if (full) dipoleDerivs[ZI][paramNum] =
						ParamDerivative.MNDODipoleDeriv((SolutionR) s,
								densityDerivs[ZI][paramNum],
								s.getRm().mats[ZI],
								paramNum);
			}
		}
		else {
			assert sPrime != null;
			HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
			if (full) dipoleDerivs[ZI][paramNum] =
					(sPrime.dipole - s.dipole) / Utils.LAMBDA;
		}
		totalGradients[ZI][paramNum] +=
				2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
		if (full) totalGradients[ZI][paramNum] +=
				800 * (s.dipole - datum[1]) * dipoleDerivs[ZI][paramNum];
	}

	@Override
	protected void computeIEDeriv(int ZI, int paramNum, Solution sPrime) {
		if (analytical) {
			if (staticDerivs[ZI][0][paramNum] != null ||
					staticDerivs[ZI][1][paramNum] != null) {
				responseDerivs[ZI][paramNum] =
						ParamDerivative.responseMatrix((SolutionR)
								s, densityDerivs[ZI][paramNum]);
				fockDerivs[ZI][paramNum] =
						staticDerivs[ZI][1][paramNum]
								.add(responseDerivs[ZI][paramNum]);
				xComplementary[ZI][paramNum] =
						ParamDerivative.xArrayComplementary((SolutionR) s,
								fockDerivs[ZI][paramNum]);
				xForIE[ZI][paramNum] =
						ParamDerivative.xarrayForIE((SolutionR) s,
								xLimited[ZI][paramNum],
								xComplementary[ZI][paramNum]);
				coeffDerivs[ZI][paramNum] = ParamDerivative
						.HOMOCoefficientDerivativeComplementary(
								xForIE[ZI][paramNum],
								(SolutionR) s);
				IEDerivs[ZI][paramNum] =
						-ParamDerivative.MNDOIEDeriv((SolutionR) s,
								coeffDerivs[ZI][paramNum],
								fockDerivs[ZI][paramNum]);
			}
		}
		else {
			assert sPrime != null;
			IEDerivs[ZI][paramNum] =
					-(sPrime.homo - s.homo) / Utils.LAMBDA;
		}
		totalGradients[ZI][paramNum] +=
				200 * -(s.homo + datum[2]) * IEDerivs[ZI][paramNum];
	}

	@Override
	protected Solution constructSPrime(int ZI, int paramNum) {
		return new SolutionR(
				Utils.perturbAtomParams(s.atoms, s.getRm().mats[ZI],
						paramNum), s.charge);
	}

	@Override
	protected SolutionR constructSExpPrime(int ZI, int paramNum) {
		return new SolutionR(
				Utils.perturbAtomParams(sExp.atoms, sExp.getRm().mats[ZI],
						paramNum), sExp.charge);
	}

	@Override
	protected double findGrad(Solution sExpPrime, int i, int xyz) {
		return GeometryDerivative.grad((SolutionR) sExpPrime, i, xyz);
	}
}
