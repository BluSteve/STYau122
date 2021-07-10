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

    public void computeHessians() {
        hessian = new double[s.getUniqueZs().length * NDDOSolution.maxParamNum]
                [s.getUniqueZs().length * NDDOSolution.maxParamNum];
        computeHessians(kind);
    }

    protected abstract void computeHessians(String kind);

    public double[][] getHessian() {
        return hessian;
    }
}
