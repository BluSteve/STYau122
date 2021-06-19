package nddoparam.param;

import nddoparam.*;
import org.jblas.DoubleMatrix;

// This should be the only class that calls ParamGradient
public abstract class ParamGradientAnalytical implements ErrorGettable { // TODO only works for restricted rn
    protected NDDOSolution s;
    protected String kind;
    protected int charge;
    protected ParamErrorFunction e;
    protected double[][] HFDerivs, dipoleDerivs, IEDerivs;
    protected DoubleMatrix[][] densityDerivs, xLimited, xComplementary, xForIE, coeffDerivs, responseDerivs, fockDerivs;
    protected DoubleMatrix[][][] staticDerivs;

    // kind = "hf_only" or "limited" or "complementary"
    // The purpose of ParamGradientHandler is to return a 1/2/3 x 5/8/13 x 1 ArrayList<DoubleMatrix>
    public ParamGradientAnalytical(NDDOSolution s, String kind) {
        this.s = s;
        this.kind = kind;
        this.HFDerivs = new double[NDDOSolution.maxParamNum][200];
        this.dipoleDerivs = new double[NDDOSolution.maxParamNum][200];
        this.IEDerivs = new double[NDDOSolution.maxParamNum][200];
        this.densityDerivs = new DoubleMatrix[NDDOSolution.maxParamNum][200];
        this.staticDerivs = new DoubleMatrix[NDDOSolution.maxParamNum][2][200];
        this.xLimited = new DoubleMatrix[NDDOSolution.maxParamNum][200];

        switch (kind) {
            case "hf_only":
                computeHFDerivs();
                break;
            case "limited":
                computeLimitedDerivs();
                break;
            case "complementary":
                computeComplementaryDerivs();
                break;
        }
    }


    protected abstract void computeDipoleDeriv(int Z, int paramNum);

    protected abstract void computeIEDeriv(int Z, int paramNum);

    protected abstract void computeBatchedDerivs(int Z);

    protected abstract void computeHFDerivs();

    protected abstract void computeLimitedDerivs();

    protected abstract void computeComplementaryDerivs();

    public ParamErrorFunction getE() {
        return this.e;
    }

    public abstract void constructErrors(double refHeat);
}
