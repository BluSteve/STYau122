package nddoparam.param;

import nddoparam.NDDOSolutionRestricted;
import scf.Utils;

public class ParamHessianRestricted extends ParamHessianAnalytical {
    public ParamHessianRestricted(NDDOSolutionRestricted s, String kind, double[] datum, NDDOSolutionRestricted sExp, boolean analytical) {
        super(s, kind, datum, sExp, analytical);
        g = new ParamGradientRestricted(s, kind, datum, sExp, true);
        g.setAnalytical(this.analytical);
        g.computeGradient();
    }

    public ParamHessianRestricted(ParamGradientRestricted g) {
        super(g.s, g.kind, g.datum, g.sExp, g.analytical);
        this.g = g;
    }

    @Override
    protected void constructGPrime(int ZIndex, int paramNum) {
        gPrime = new ParamGradientRestricted(new NDDOSolutionRestricted(Utils.perturbAtomParams(s.atoms,
                s.getUniqueZs()[ZIndex], paramNum), s.charge), kind, datum, (NDDOSolutionRestricted) sExp, true);
        gPrime.setAnalytical(analytical);
    }
}
