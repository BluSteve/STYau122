package nddoparam.param;

import nddoparam.NDDOSolutionRestricted;
import scf.Utils;

public class ParamHessianRestricted extends ParamHessianAnalytical {
    public ParamHessianRestricted(NDDOSolutionRestricted s, String kind, double[] datum, NDDOSolutionRestricted sExp) {
        super(s, kind, datum, sExp);
        g = new ParamGradientRestricted(s, kind, datum, sExp);
        g.setAnalytical(analytical);
        g.computeDerivs();
    }


    @Override
    protected void constructGPrime(int ZIndex, int paramNum) {
        gPrime = new ParamGradientRestricted(new NDDOSolutionRestricted(Utils.perturbAtomParams(s.atoms,
                s.getUniqueZs()[ZIndex], paramNum), s.charge), kind, datum, (NDDOSolutionRestricted) sExp);
    }
}
