package nddo.param;

import nddo.Constants;
import nddo.geometry.GeometryDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionU;
import tools.Utils;

class ParamGradientU extends ParamGradient {
	private static final String errorMessage =
			"Analytical derivatives have yet to be implemented for " +
					"unrestricted!";

	protected ParamGradientU(SolutionU s, double[] datum, SolutionU sExp,
							 boolean analytical) {
		super(s, datum, sExp, analytical);
	}

	@Override
	protected void computeBatchedDerivs(int firstZIndex, int firstParamNum) {
		s.rm.getLogger().error(errorMessage);
	}

	@Override
	protected void computeHFDeriv(int ZI, int paramNum, Solution sPrime) {
		if (analytical) {
			s.rm.getLogger().error(errorMessage);
		}
		else {
			HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Constants.LAMBDA;
		}

		totalGradients[ZI][paramNum] +=
				2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
	}

	@Override
	protected void computeDipoleDeriv(int ZI, int paramNum, boolean full,
									  Solution sPrime) {
		if (analytical) {
			s.rm.getLogger().error(errorMessage);
		}
		else {
			HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Constants.LAMBDA;
			if (full) dipoleDerivs[ZI][paramNum] =
					(sPrime.dipole - s.dipole) / Constants.LAMBDA;
		}
		if (full) totalGradients[ZI][paramNum] +=
				800 * (s.dipole - datum[1]) * dipoleDerivs[ZI][paramNum];
	}

	@Override
	protected void computeIEDeriv(int ZI, int paramNum, Solution sPrime) {
		if (analytical) {
			s.rm.getLogger().error(errorMessage);
		}
		else {
			IEDerivs[ZI][paramNum] = -(sPrime.homo - s.homo) / Constants.LAMBDA;
		}
		totalGradients[ZI][paramNum] +=
				200 * -(s.homo + datum[2]) * IEDerivs[ZI][paramNum];
	}

	@Override
	protected Solution constructSPrime(int ZI, int paramNum) {
		return s.withNewAtoms(
				Utils.perturbAtomParams(s.atoms, s.rm.mats[ZI],
						paramNum));
	}

	@Override
	protected Solution constructSExpPrime(int ZI, int paramNum) {
		return sExp.withNewAtoms(
				Utils.perturbAtomParams(sExp.atoms, sExp.rm.mats[ZI],
						paramNum));
	}

	@Override
	protected double findGrad(Solution sExpPrime, int i, int xyz) {
		return GeometryDerivative.grad((SolutionU) sExpPrime, i, xyz);
	}
}
