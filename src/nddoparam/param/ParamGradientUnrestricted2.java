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
    protected void constructSPrime(int Z, int paramNum) {
        sPrime = new NDDOSolutionUnrestricted(Utils.perturbAtomParams(s.atoms, s.getUniqueZs()[Z], paramNum), s.charge,
                s.multiplicity);
    }

    // Compiles all necessary fock matrices into one array before using the Pople algorithm, for faster computation.
    // This function is the only thing that's not computed on a Z, paramNum level.
    // Will not compute anything before the firstZIndex and the firstParamIndex.
    @Override
    protected void computeBatchedDerivs(int firstZIndex, int firstParamIndex) {

    }

    @Override
    protected void computeHFDeriv(int Z, int paramNum) {
        if (analytical) {
        } else
            HFDerivs[Z][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
        totalGradients[Z][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[Z][paramNum];
    }

    // `full` sets whether dipole itself is actually computed
    // also computes HF
    @Override
    protected void computeDipoleDeriv(int Z, int paramNum, boolean full) {
        if (analytical) {

        } else {
            HFDerivs[Z][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
            if (full) dipoleDerivs[Z][paramNum] = (sPrime.dipole - s.dipole) / Utils.LAMBDA;
        }
        totalGradients[Z][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[Z][paramNum];
        if (full) totalGradients[Z][paramNum] += 800 * (s.dipole - datum[1]) * dipoleDerivs[Z][paramNum];
    }

    @Override
    protected void computeIEDeriv(int Z, int paramNum) {
        if (analytical) {

        } else {
            IEDerivs[Z][paramNum] = -(sPrime.homo - s.homo) / Utils.LAMBDA;
        }
        totalGradients[Z][paramNum] += 200 * (s.homo + datum[2]) * IEDerivs[Z][paramNum];
    }

    @Override
    protected void computeGeomDeriv(int Z, int paramNum) {
        sExpPrime = new NDDOSolutionUnrestricted(Utils.perturbAtomParams(sExp.atoms, s.getUniqueZs()[Z], paramNum), s.charge,
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
        geomDerivs[Z][paramNum] = 1 / LAMBDA * (geomGradient - e.geomGradient);
        totalGradients[Z][paramNum] += 0.000049 * e.geomGradient * geomDerivs[Z][paramNum];
    }
}
