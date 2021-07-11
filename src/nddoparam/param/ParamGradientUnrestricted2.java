package nddoparam.param;

import nddoparam.NDDODerivative;
import nddoparam.NDDOSolutionUnrestricted;
import scf.Utils;

public class ParamGradientUnrestricted2 extends ParamGradientAnalytical {

    public ParamGradientUnrestricted2(NDDOSolutionUnrestricted s, String kind, double[] datum,
                                      NDDOSolutionUnrestricted sExp, boolean analytical) {
        super(s, kind, datum, sExp, analytical);
        e = new ParamErrorFunctionUnrestricted(s, datum[0]);
        if (datum[1] != 0) e.addDipoleError(datum[1]);
        if (datum[2] != 0) e.addIEError(datum[2]);

        if (this.sExp != null) {
            isExpAvail = true;
            e.createExpGeom(this.sExp.atoms, this.sExp);
            e.addGeomError();
        }
    }

    @Override
    protected void constructSPrime(int ZI, int paramNum) {
        sPrime = new NDDOSolutionUnrestricted(Utils.perturbAtomParams(s.atoms, s.getUniqueZs()[ZI], paramNum), s.charge,
                s.multiplicity);
    }

    @Override
    protected void computeBatchedDerivs(int firstZIndex, int firstParamIndex) {

    }

    @Override
    protected void computeHFDeriv(int ZI, int paramNum) {
        if (analytical) {
        } else
            HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
        totalGradients[ZI][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
    }

    @Override
    protected void computeDipoleDeriv(int ZI, int paramNum, boolean full) {
        if (analytical) {

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

        } else {
            IEDerivs[ZI][paramNum] = -(sPrime.homo - s.homo) / Utils.LAMBDA;
        }
        totalGradients[ZI][paramNum] += 200 * (s.homo + datum[2]) * IEDerivs[ZI][paramNum];
    }

    @Override
    protected void computeGeomDeriv(int ZI, int paramNum) {
        sExpPrime = new NDDOSolutionUnrestricted(Utils.perturbAtomParams(sExp.atoms, s.getUniqueZs()[ZI], paramNum), s.charge,
                s.multiplicity);
        double sum = 0;
        double d;
        for (int i = 0; i < sExpPrime.atoms.length; i++) {
            for (int j = 0; j < 3; j++) {
                d = NDDODerivative.grad((NDDOSolutionUnrestricted) sExpPrime, i, j);
                sum += d * d;
            }
        }
        double geomGradient = 627.5 * Math.sqrt(sum);
        geomDerivs[ZI][paramNum] = 1 / LAMBDA * (geomGradient - e.geomGradient);
        totalGradients[ZI][paramNum] += 0.000049 * e.geomGradient * geomDerivs[ZI][paramNum];
    }
}
