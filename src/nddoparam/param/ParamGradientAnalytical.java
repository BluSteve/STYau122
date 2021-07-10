package nddoparam.param;

import nddoparam.NDDOSolution;
import org.jblas.DoubleMatrix;
import scf.Utils;

public abstract class ParamGradientAnalytical implements ErrorGettable { // TODO only works for restricted rn
    protected NDDOSolution s, sPrime, sExp;
    protected ParamErrorFunction e;
    protected String kind;
    protected boolean isExpAvail;
    protected double[] datum;
    protected double[][] HFDerivs, dipoleDerivs, IEDerivs, geomDerivs, totalDerivs;
    protected DoubleMatrix[][] densityDerivs, xLimited, xComplementary, xForIE, coeffDerivs, responseDerivs, fockDerivs;
    protected DoubleMatrix[][][] staticDerivs;
    protected static final double LAMBDA = 1E-7;

    // kind = "hf_only" or "limited" or "complementary"
    public ParamGradientAnalytical(NDDOSolution s, String kind, double[] datum, NDDOSolution sExp) {
        this.s = s;
        this.kind = kind;
        this.datum = datum;
        this.sExp = sExp;
    }

    public void computeDerivs() {
        geomDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
        totalDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
        switch (kind) {
            case "a":
                HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                computeADerivs();
                break;
            case "b":
                HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                dipoleDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                densityDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                staticDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum][2];
                xLimited = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                computeBDerivs();
                break;
            case "c":
                HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                IEDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                densityDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                staticDerivs = new DoubleMatrix[s.getUniqueZs().length][2][NDDOSolution.maxParamNum];
                xLimited = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                xComplementary = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                xForIE = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                coeffDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                responseDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                fockDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                computeCDerivs();
                break;
            case "d":
                HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                dipoleDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                IEDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                densityDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                staticDerivs = new DoubleMatrix[s.getUniqueZs().length][2][NDDOSolution.maxParamNum];
                xLimited = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                xComplementary = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                xForIE = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                coeffDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                responseDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                fockDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                computeDDerivs();
                break;
        }
    }

    protected abstract void computeGeomDeriv(int Z, int paramNum);

    protected abstract void computeHFDeriv(int Z, int paramNum);

    protected abstract void computeDipoleDeriv(int Z, int paramNum, boolean lite);

    protected abstract void computeIEDeriv(int Z, int paramNum);

    protected abstract void computeBatchedDerivs(int firstZIndex, int firstParamIndex);

    protected abstract void computeADerivs();

    protected abstract void computeBDerivs();

    protected abstract void computeCDerivs();

    protected abstract void computeDDerivs();

    public ParamErrorFunction getE() {
        return this.e;
    }

    public double[][] getHFDerivs() {
        return HFDerivs;
    }

    public double[][] getDipoleDerivs() {
        return dipoleDerivs;
    }

    public double[][] getIEDerivs() {
        return IEDerivs;
    }

    public double[][] getGeomDerivs() {
        return geomDerivs;
    }

    public double[][] getTotalDerivs() {
        return totalDerivs;
    }

    public NDDOSolution getS() {
        return s;
    }
}
