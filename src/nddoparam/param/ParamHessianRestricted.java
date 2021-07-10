package nddoparam.param;

import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;

public class ParamHessianRestricted extends ParamHessianAnalytical {
    public ParamHessianRestricted(NDDOSolutionRestricted s, String kind, int charge, double[] datum, NDDOSolution sExp) {
        super(s, kind, charge, datum, sExp);
        g = new ParamGradientRestricted(s, kind, charge, datum, sExp);
        g.computeDerivs();
    }


    @Override
    public ParamErrorFunction getE() {
        return null;
    }

    @Override
    protected void computeAHessians() {
        g.getTotalDerivs();
    }

    @Override
    protected void computeBHessians() {

    }

    @Override
    protected void computeCHessians() {

    }

    @Override
    protected void computeDHessians() {

    }
}
