package nddoparam.param;

import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;
import scf.Utils;

public class ParamHessianRestricted extends ParamHessianAnalytical {
    public ParamHessianRestricted(NDDOSolution s, String kind, double[] datum, NDDOSolution sExp) {
        super(s, kind, datum, sExp);
        g = new ParamGradientRestricted(s, kind, datum, sExp);
        g.computeDerivs();
    }


    @Override
    public ParamErrorFunction getE() {
        return g.getE();
    }

    @Override
    protected void computeHessians(String kind) {
        for (int ZIndex2 = 0; ZIndex2 < s.getUniqueZs().length; ZIndex2++) {
            for (int paramNum2 = 0; paramNum2 < NDDOSolution.maxParamNum; paramNum2++) {
                gPrime = new ParamGradientRestricted(new NDDOSolutionRestricted(Utils.perturbAtomParams(s.atoms, paramNum2,
                        s.getUniqueZs()[ZIndex2]), s.charge), kind, datum, sExp);
                if (kind.equals("b") || kind.equals("c") || kind.equals("d"))
                    gPrime.computeBatchedDerivs(ZIndex2, paramNum2);

                for (int ZIndex1 = ZIndex2; ZIndex1 < s.getUniqueZs().length; ZIndex1++) {
                    for (int paramNum1 = paramNum2; paramNum1 < NDDOSolution.maxParamNum; paramNum1++) {
                        switch (kind) {
                            case "a":
                                gPrime.computeHFDeriv(ZIndex1, paramNum1);
                                break;
                            case "b":
                                gPrime.computeDipoleDeriv(ZIndex1, paramNum1, true);
                                break;
                            case "c":
                                gPrime.computeDipoleDeriv(ZIndex1, paramNum1, false);
                                gPrime.computeIEDeriv(ZIndex1, paramNum1);
                                break;
                            case "d":
                                gPrime.computeDipoleDeriv(ZIndex1, paramNum1, true);
                                gPrime.computeIEDeriv(ZIndex1, paramNum1);
                                break;
                        }
                        if (gPrime.isExpAvail) gPrime.computeGeomDeriv(ZIndex1, paramNum1);
                        hessian[ZIndex2 * NDDOSolution.maxParamNum + paramNum2]
                                [ZIndex1 * NDDOSolution.maxParamNum + paramNum1] = (gPrime.getGradients()[ZIndex1][paramNum1]
                                - g.getGradients()[ZIndex1][paramNum1]) / Utils.lambda;
                        hessian[ZIndex1 * NDDOSolution.maxParamNum + paramNum1]
                                [ZIndex2 * NDDOSolution.maxParamNum + paramNum2] =
                                hessian[ZIndex2 * NDDOSolution.maxParamNum + paramNum2]
                                        [ZIndex1 * NDDOSolution.maxParamNum + paramNum1];
                    }
                }
            }
        }
    }
}
