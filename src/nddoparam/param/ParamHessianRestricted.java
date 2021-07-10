package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;
import scf.Utils;

import java.util.Arrays;

public class ParamHessianRestricted extends ParamHessianAnalytical {
    public ParamHessianRestricted(NDDOSolutionRestricted s, String kind, double[] datum, NDDOSolution sExp) {
        super(s, kind, datum, sExp);
        g = new ParamGradientRestricted(s, kind, datum, sExp);
        g.computeDerivs();
    }


    @Override
    public ParamErrorFunction getE() {
        return null;
    }

    @Override
    protected void computeAHessians() {
    }

    @Override
    protected void computeAHessian(int ZIndex2, int paramNum2) {
        computeBatchedDerivs(ZIndex2, paramNum2);
        for (int ZIndex1 = ZIndex2; ZIndex1 < s.getUniqueZs().length; ZIndex1++) {
            for (int paramNum1 = paramNum2; paramNum1 < NDDOSolution.maxParamNum; paramNum1++) {
                gPrime.computeHFDeriv(s.getUniqueZs()[ZIndex1], paramNum1);
                if (gPrime.isExpAvail) gPrime.computeGeomDeriv(s.getUniqueZs()[ZIndex1], paramNum1);
                hessian[ZIndex2 * NDDOSolution.maxParamNum + paramNum2]
                        [ZIndex1 * NDDOSolution.maxParamNum + paramNum1] = gPrime.getTotalDerivs()[s.getUniqueZs()[ZIndex1]][paramNum1]
                        - g.getTotalDerivs()[s.getUniqueZs()[ZIndex1]][paramNum1];
            }
        }
    }

    @Override
    protected void computeBatchedDerivs(int ZIndex, int paramNum2) {
        gPrime = new ParamGradientRestricted(new NDDOSolutionRestricted(Utils.perturbAtomParams(s.atoms, paramNum2,
                s.getUniqueZs()[ZIndex]), s.charge), kind, datum, sExp);
        gPrime.computeBatchedDerivs(s.getUniqueZs()[ZIndex], paramNum2);
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
