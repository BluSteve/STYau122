package nddoparam.param;

import nddoparam.GeometryDerivative;
import nddoparam.Solution;
import nddoparam.SolutionU;
import scf.Utils;

public class ParamGradientU extends ParamGradient {
	private final String errorMessage =
			"Analytical derivatives have yet to be implemented for " +
					"unrestricted!";

	public ParamGradientU(SolutionU s, double[] datum, SolutionU sExp,
						  boolean analytical) {
		super(s, datum, sExp, analytical);
		initializeArrays();

		e = new ParamErrorFunctionU(s, datum[0]);
		errorFunctionRoutine();
	}

	@Override
	protected void computeBatchedDerivs(int firstZIndex, int firstParamIndex) {
		System.err.println(errorMessage);
	}

	@Override
	protected void computeHFDeriv(int ZI, int paramNum,
								  Solution sPrime) {
		if (analytical) {
			System.err.println(errorMessage);
		}
		else {
			HFDerivs[ZI][paramNum] = (this.sPrime.hf - s.hf) / Utils.LAMBDA;
		}

		totalGradients[ZI][paramNum] +=
				2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
	}

	@Override
	protected void computeDipoleDeriv(int ZI, int paramNum, boolean full,
									  Solution sPrime) {
		if (analytical) {
			System.err.println(errorMessage);
		}
		else {
			HFDerivs[ZI][paramNum] = (this.sPrime.hf - s.hf) / Utils.LAMBDA;
			if (full) dipoleDerivs[ZI][paramNum] =
					(this.sPrime.dipole - s.dipole) / Utils.LAMBDA;
		}
		totalGradients[ZI][paramNum] +=
				2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
		if (full) totalGradients[ZI][paramNum] +=
				800 * (s.dipole - datum[1]) * dipoleDerivs[ZI][paramNum];
	}

	@Override
	protected void computeIEDeriv(int ZI, int paramNum,
								  Solution sPrime) {
		if (analytical) {
			System.err.println(errorMessage);
		}
		else {
			IEDerivs[ZI][paramNum] = -(this.sPrime.homo - s.homo) / Utils.LAMBDA;
		}
		totalGradients[ZI][paramNum] +=
				200 * -(s.homo + datum[2]) * IEDerivs[ZI][paramNum];
	}

	@Override
	protected Solution constructSPrime(int ZI, int paramNum) {
		return new SolutionU(
				Utils.perturbAtomParams(s.atoms, s.getUniqueZs()[ZI],
						paramNum), s.charge, s.multiplicity);
	}

	@Override
	protected SolutionU constructSExpPrime(int Z, int paramNum) {
		SolutionU sExpPrime = new SolutionU(
				Utils.perturbAtomParams(sExp.atoms, sExp.getUniqueZs()[Z],
						paramNum),
				sExp.charge,
				sExp.multiplicity);
		return sExpPrime;
	}

	@Override
	protected double findGrad(Solution sExpPrime, int i, int j) {
		return GeometryDerivative.grad((SolutionU) sExpPrime, i, j);
	}
}
