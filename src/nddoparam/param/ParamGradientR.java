package nddoparam.param;

import nddoparam.GeometryDerivative;
import nddoparam.Solution;
import nddoparam.SolutionR;
import org.jblas.DoubleMatrix;
import scf.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class ParamGradientR extends ParamGradient {

	public ParamGradientR(SolutionR s, double[] datum, SolutionR sExp,
						  boolean analytical) {
		super(s, datum, sExp, analytical);
		initializeArrays();

		e = new ParamErrorFunctionR(s, datum[0]);
		errorFunctionRoutine();
	}

	// Compiles all necessary fock matrices into one array before using the
	// Pople algorithm, for faster computation. This function is the only thing
	// that's not computed on a Z, paramNum level. Will not compute anything
	// before the firstZIndex and the firstParamIndex.
	@Override
	protected void computeBatchedDerivs(int firstZIndex, int firstParamIndex) {
		ArrayList<DoubleMatrix> aggregate =
				new ArrayList<>(s.getUniqueZs().length * Solution.maxParamNum);
		for (int Z = firstZIndex; Z < s.getUniqueZs().length; Z++) {
			if (Z == firstZIndex)
				staticDerivs[Z] = ParamDerivative
						.MNDOStaticMatrixDeriv((SolutionR) s,
								s.getUniqueZs()[Z],
								firstParamIndex);
			else staticDerivs[Z] = ParamDerivative
					.MNDOStaticMatrixDeriv((SolutionR) s,
							s.getUniqueZs()[Z],
							0);
			Collections.addAll(aggregate, staticDerivs[Z][1]);
		}

		DoubleMatrix[] aggregateArray = new DoubleMatrix[aggregate.size()];
		for (int x = 0; x < aggregate.size(); x++) {
			aggregateArray[x] = aggregate.get(x);
		}
		DoubleMatrix[] xLimitedAggregate =
				ParamDerivative.xArrayLimitedPople((SolutionR) s,
						aggregateArray);
		int i = 0;
		for (int Z = firstZIndex; Z < s.getUniqueZs().length; Z++) {
			xLimited[Z] =
					Arrays.copyOfRange(xLimitedAggregate,
							i * Solution.maxParamNum,
							i * Solution.maxParamNum + Solution.maxParamNum);
			i++;
		}
	}

	@Override
	protected void computeHFDeriv(int ZI, int paramNum) {
		if (analytical)
			HFDerivs[ZI][paramNum] = ParamDerivative
					.HFDeriv((SolutionR) s, s.getUniqueZs()[ZI], paramNum);
		else
			HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
		totalGradients[ZI][paramNum] +=
				2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
	}

	// `full` sets whether dipole itself is actually computed
	// also computes HF
	@Override
	protected void computeDipoleDeriv(int ZI, int paramNum, boolean full) {
		if (analytical) {
			// alpha and eisol have to be computed apart from the rest of
			// HFDerivs. I.e not as a part of dipole.
			if (paramNum == 0 || paramNum == 7) {
				HFDerivs[ZI][paramNum] = ParamDerivative
						.HFDeriv((SolutionR) s, s.getUniqueZs()[ZI],
								paramNum);
			}
			if (staticDerivs[ZI][0][paramNum] != null ||
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
								s.getUniqueZs()[ZI],
								paramNum);
			}
		}
		else {
			HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
			if (full)
				dipoleDerivs[ZI][paramNum] =
						(sPrime.dipole - s.dipole) / Utils.LAMBDA;
		}
		totalGradients[ZI][paramNum] +=
				2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
		if (full) totalGradients[ZI][paramNum] +=
				800 * (s.dipole - datum[1]) * dipoleDerivs[ZI][paramNum];
	}

	@Override
	protected void computeIEDeriv(int ZI, int paramNum) {
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
			IEDerivs[ZI][paramNum] = -(sPrime.homo - s.homo) / Utils.LAMBDA;
		}
		totalGradients[ZI][paramNum] +=
				200 * -(s.homo + datum[2]) * IEDerivs[ZI][paramNum];
	}

	@Override
	protected void constructSPrime(int ZI, int paramNum) {
		sPrime = new SolutionR(
				Utils.perturbAtomParams(s.atoms, s.getUniqueZs()[ZI],
						paramNum),
				s.charge);
	}

	@Override
	protected void constructSExpPrime(int Z, int paramNum) {
		sExpPrime = new SolutionR(
				Utils.perturbAtomParams(sExp.atoms, s.getUniqueZs()[Z],
						paramNum),
				s.charge);
	}

	@Override
	protected double findGrad(int i, int j) {
		return GeometryDerivative.grad((SolutionR) sExpPrime, i, j);
	}
}