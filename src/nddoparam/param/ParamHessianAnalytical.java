package nddoparam.param;

import nddoparam.NDDOSolution;
import scf.Utils;

public abstract class ParamHessianAnalytical implements ErrorGettable {
    protected ParamGradientAnalytical g, gPrime;
    protected NDDOSolution s, sExp;
    protected String kind;
    protected double[] datum;
    protected double[][] hessian;

    public ParamHessianAnalytical(NDDOSolution s, String kind, double[] datum, NDDOSolution sExp) {
        this.s = s;
        this.kind = kind;
        this.datum = datum;
        this.sExp = sExp;
        hessian = new double[s.getUniqueZs().length * NDDOSolution.maxParamNum]
                [s.getUniqueZs().length * NDDOSolution.maxParamNum];
    }

    protected void computeHessian() {
        for (int ZIndex2 = 0; ZIndex2 < s.getUniqueZs().length; ZIndex2++) {
            for (int paramNum2 : s.getNeededParams()[s.getUniqueZs()[ZIndex2]]) {
                constructGPrime(ZIndex2, paramNum2);
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

    protected abstract void constructGPrime(int ZIndex, int paramNum);

    public double[][] getHessian() {
        return hessian;
    }

    @Override
    public ParamErrorFunction getE() {
        return g.getE();
    }
}
