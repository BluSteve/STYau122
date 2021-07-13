package nddoparam.param;

import nddoparam.NDDODerivative;
import nddoparam.NDDOSolutionUnrestricted;
import scf.Utils;

public class ParamGradientUnrestricted2 extends ParamGradientAnalytical {
    private final String errorMessage = "Analytical derivatives have yet to be implemented for unrestricted!";

    public ParamGradientUnrestricted2(NDDOSolutionUnrestricted s, String kind, double[] datum,
                                      NDDOSolutionUnrestricted sExp, boolean analytical) {
        super(s, kind, datum, sExp, analytical);
        initializeArrays();

        e = new ParamErrorFunctionUnrestricted(s, datum[0]);
        errorFunctionRoutine();
    }

    @Override
    protected void computeBatchedDerivs(int firstZIndex, int firstParamIndex) {
        System.err.println(errorMessage);
    }

    @Override
    protected void computeHFDeriv(int ZI, int paramNum) {
        if (analytical) {
            System.err.println(errorMessage);
        } else
            HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
        totalGradients[ZI][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
    }

    @Override
    protected void computeDipoleDeriv(int ZI, int paramNum, boolean full) {
        if (analytical) {
            System.err.println(errorMessage);
        } else {
            HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
            if (full) dipoleDerivs[ZI][paramNum] = (sPrime.dipole - s.dipole) / Utils.LAMBDA;
        }
        totalGradients[ZI][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
        if (full) totalGradients[ZI][paramNum] += 800 * (s.dipole - datum[1]) * dipoleDerivs[ZI][paramNum];
    }

    @Override
    protected void computeIEDeriv(int ZI, int paramNum) {
        if (analytical) {
            System.err.println(errorMessage);
        } else {
            IEDerivs[ZI][paramNum] = -(sPrime.homo - s.homo) / Utils.LAMBDA;
        }
        totalGradients[ZI][paramNum] += 200 * -(s.homo + datum[2]) * IEDerivs[ZI][paramNum];
    }

    @Override
    protected void constructSPrime(int ZI, int paramNum) {
        sPrime = new NDDOSolutionUnrestricted(Utils.perturbAtomParams(s.atoms, s.getUniqueZs()[ZI], paramNum), s.charge,
                s.multiplicity);
    }

    @Override
    protected void constructSExpPrime(int Z, int paramNum) {
        // TODO sExp.charge or s.charge?
        sExpPrime = new NDDOSolutionUnrestricted(Utils.perturbAtomParams(sExp.atoms, s.getUniqueZs()[Z], paramNum), s.charge,
                s.multiplicity);
    }

    @Override
    protected double findGrad(int i, int j) {
        return NDDODerivative.grad((NDDOSolutionUnrestricted) sExpPrime, i, j);
    }
}
