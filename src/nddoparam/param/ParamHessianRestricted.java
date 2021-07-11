package nddoparam.param;

import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;
import scf.Utils;

public class ParamHessianRestricted extends ParamHessianAnalytical {
    public ParamHessianRestricted(NDDOSolution s, String kind, double[] datum, NDDOSolution sExp) {
        super(s, kind, datum, sExp);
        g = new ParamGradientRestricted(s, kind, datum, sExp, true);
        g.computeDerivs();
    }


    @Override
    protected void computeHessian(String kind) {
        for (int ZIndex2 = 0; ZIndex2 < s.getUniqueZs().length; ZIndex2++) {
            for (int paramNum2: s.getNeededParams()[s.getUniqueZs()[ZIndex2]]) {
                gPrime = new ParamGradientRestricted(new NDDOSolutionRestricted(Utils.perturbAtomParams(s.atoms,
                        s.getUniqueZs()[ZIndex2], paramNum2), s.charge), kind, datum, sExp, true);
                if (kind.equals("b") || kind.equals("c") || kind.equals("d"))
                    gPrime.computeBatchedDerivs(ZIndex2, paramNum2);


                boolean needed;
                for (int ZIndex1 = ZIndex2; ZIndex1 < s.getUniqueZs().length; ZIndex1++) {
                    for (int paramNum1 = paramNum2; paramNum1 < NDDOSolution.maxParamNum; paramNum1++) {
                        needed = false;
                        for (int p : s.getNeededParams()[s.getUniqueZs()[ZIndex1]]) {
                            if (paramNum1 == p) {
                                needed = true;
                                break;
                            }
                        }


                        if (needed) {
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
                                    [ZIndex1 * NDDOSolution.maxParamNum + paramNum1] =
                                    (gPrime.getTotalGradients()[ZIndex1][paramNum1]
                                    - g.getTotalGradients()[ZIndex1][paramNum1]) / Utils.LAMBDA;
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
}
