package nddoparam.param;

import nddoparam.NDDOSolutionUnrestricted;
import scf.Utils;

public class ParamHessianUnrestricted2 extends ParamHessianAnalytical {
    public ParamHessianUnrestricted2(NDDOSolutionUnrestricted s, String kind, double[] datum, NDDOSolutionUnrestricted sExp, boolean analytical) {
        // TODO Change this line's hardcoding once analytical has been implemented for unrestricted.
        super(s, kind, datum, sExp, false);
        g = new ParamGradientUnrestricted2(s, kind, datum, sExp, analytical);
        g.computeGradients();
    }

    public ParamHessianUnrestricted2(ParamGradientUnrestricted2 g, boolean analytical) {
        // TODO Change this line's hardcoding once analytical has been implemented for unrestricted.
        super(g.s, g.kind, g.datum, g.sExp, false);
        this.g = g;
    }

    @Override
    protected void constructGPrime(int ZIndex, int paramNum) {
        gPrime = new ParamGradientUnrestricted2(new NDDOSolutionUnrestricted(Utils.perturbAtomParams(s.atoms,
                s.getUniqueZs()[ZIndex], paramNum), s.charge, s.multiplicity), kind, datum, (NDDOSolutionUnrestricted) sExp, analytical);
    }
}
