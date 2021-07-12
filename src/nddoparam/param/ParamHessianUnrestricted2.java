package nddoparam.param;

import nddoparam.NDDOSolutionUnrestricted;
import scf.Utils;

public class ParamHessianUnrestricted2 extends ParamHessianAnalytical {
    public ParamHessianUnrestricted2(NDDOSolutionUnrestricted s, String kind, double[] datum, NDDOSolutionUnrestricted sExp, boolean analytical) {
        super(s, kind, datum, sExp, analytical);
        g = new ParamGradientUnrestricted2(s, kind, datum, sExp, true);
        g.setAnalytical(this.analytical);
        g.computeGradient();
    }

    public ParamHessianUnrestricted2(ParamGradientUnrestricted2 g) {
        super(g.s, g.kind, g.datum, g.sExp, g.analytical);
        this.g = g;
    }

    @Override
    protected void constructGPrime(int ZIndex, int paramNum) {
        gPrime = new ParamGradientUnrestricted2(new NDDOSolutionUnrestricted(Utils.perturbAtomParams(s.atoms,
                s.getUniqueZs()[ZIndex], paramNum), s.charge, s.multiplicity), kind, datum, (NDDOSolutionUnrestricted) sExp, true);
        gPrime.setAnalytical(analytical);
    }
}
