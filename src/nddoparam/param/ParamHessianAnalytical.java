package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;
import org.jblas.DoubleMatrix;
import scf.Utils;

public abstract class ParamHessianAnalytical implements ErrorGettable{
    protected ParamGradientAnalytical g;
    protected NDDOSolution s, sExp;
    protected String kind;
    protected int charge;
    protected double[] datum;

    public ParamHessianAnalytical(NDDOSolution s, String kind, int charge, double[] datum, NDDOSolution sExp) {
        this.s = s;
        this.kind = kind;
        this.charge = charge;
        this.datum = datum;
        this.sExp = sExp;
    }

    public void computeDerivs() {
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
    protected abstract void computeBHessians();
    protected abstract void computeCHessians();
    protected abstract void computeDHessians();
}
