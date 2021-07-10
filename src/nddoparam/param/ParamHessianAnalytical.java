package nddoparam.param;

import nddoparam.NDDOSolution;

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
    }

    public void computeDerivs() {
        hessian = new double[s.getUniqueZs().length * NDDOSolution.maxParamNum]
                [s.getUniqueZs().length * NDDOSolution.maxParamNum];
        switch (kind) {
            case "a":
                computeAHessians();
                break;
            case "b":
                computeBHessians();
                break;
            case "c":
                computeCHessians();
                break;
            case "d":
                computeDHessians();
                break;
        }
    }

    protected abstract void computeAHessians();

    protected abstract void computeAHessian(int ZIndex, int paramNum2);

    protected abstract void computeBatchedDerivs(int ZIndex, int paramNum2);

    protected abstract void computeBHessians();

    protected abstract void computeCHessians();

    protected abstract void computeDHessians();
}
