package nddoparam.param;

import nddoparam.NDDOSolutionRestricted;
import nddoparam.NDDOSolutionUnrestricted;
import scf.Utils;

public class ParamHessianUnrestricted2 extends ParamHessianAnalytical {
    public ParamHessianUnrestricted2(NDDOSolutionUnrestricted s, String kind, double[] datum, NDDOSolutionUnrestricted sExp) {
        super(s, kind, datum, sExp);
        g = new ParamGradientUnrestricted2(s, kind, datum, sExp);
        g.setAnalytical(analytical);
        g.computeDerivs();
    }


    @Override
    protected void constructGPrime(int ZIndex, int paramNum) {
        gPrime = new ParamGradientUnrestricted2(new NDDOSolutionUnrestricted(Utils.perturbAtomParams(s.atoms,
                s.getUniqueZs()[ZIndex], paramNum), s.charge, s.multiplicity), kind, datum, (NDDOSolutionUnrestricted) sExp);
    }
}
