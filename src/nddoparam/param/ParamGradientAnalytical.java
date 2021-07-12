package nddoparam.param;

import nddoparam.NDDOSolution;
import org.jblas.DoubleMatrix;

public abstract class ParamGradientAnalytical implements ErrorGettable {
    protected NDDOSolution s, sPrime, sExpPrime, sExp;
    protected ParamErrorFunction e;
    protected String kind;
    protected boolean isExpAvail, analytical;

    public boolean isAnalytical() {
        return analytical;
    }

    public void setAnalytical(boolean analytical) {
        this.analytical = analytical;
    }

    protected double[] datum;
    protected double[][] HFDerivs, dipoleDerivs, IEDerivs, geomDerivs, totalGradients;
    protected DoubleMatrix[][] densityDerivs, xLimited, xComplementary, xForIE, coeffDerivs, responseDerivs, fockDerivs;
    protected DoubleMatrix[][][] staticDerivs;
    protected static final double LAMBDA = 1E-7;

    public ParamGradientAnalytical(NDDOSolution s, String kind, double[] datum, NDDOSolution sExp) {
        this.s = s;
        this.kind = kind;
        this.datum = datum;
        this.sExp = sExp;
        this.analytical = true;

        geomDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
        totalGradients = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
        switch (kind) {
            case "a":
                HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                break;
            case "b":
                HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                dipoleDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                if (this.analytical) {
                    densityDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    staticDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum][2];
                    xLimited = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                }
                break;
            case "c":
                HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                IEDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                if (this.analytical) {
                    densityDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    staticDerivs = new DoubleMatrix[s.getUniqueZs().length][2][NDDOSolution.maxParamNum];
                    xLimited = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                    xComplementary = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    xForIE = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    coeffDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    responseDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    fockDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                }
                break;
            case "d":
                HFDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                dipoleDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                IEDerivs = new double[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                if (this.analytical) {
                    densityDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    staticDerivs = new DoubleMatrix[s.getUniqueZs().length][2][NDDOSolution.maxParamNum];
                    xLimited = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];

                    xComplementary = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    xForIE = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    coeffDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    responseDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                    fockDerivs = new DoubleMatrix[s.getUniqueZs().length][NDDOSolution.maxParamNum];
                }
                break;
        }
    }

    public void computeDerivs() {
        if (analytical && (kind.equals("b") || kind.equals("c") || kind.equals("d")))
            computeBatchedDerivs(0, 0);
        for (int Z = 0; Z < s.getUniqueZs().length; Z++) {
            for (int paramNum : s.getNeededParams()[s.getUniqueZs()[Z]]) {
                if (!analytical) constructSPrime(Z, paramNum);
                switch (kind) {
                    case "a":
                        computeHFDeriv(Z, paramNum);
                        break;
                    case "b":
                        computeDipoleDeriv(Z, paramNum, true);
                        break;
                    case "c":
                        computeDipoleDeriv(Z, paramNum, false);
                        computeIEDeriv(Z, paramNum);
                        break;
                    case "d":
                        computeDipoleDeriv(Z, paramNum, true);
                        computeIEDeriv(Z, paramNum);
                        break;
                }
                if (isExpAvail) computeGeomDeriv(Z, paramNum);
            }
        }
    }

    protected abstract void constructSPrime(int Z, int paramNum);

    protected abstract void computeBatchedDerivs(int firstZIndex, int firstParamIndex);

    protected abstract void computeHFDeriv(int Z, int paramNum);

    protected abstract void computeDipoleDeriv(int Z, int paramNum, boolean full);

    protected abstract void computeIEDeriv(int Z, int paramNum);

    protected abstract void computeGeomDeriv(int Z, int paramNum);

    public ParamErrorFunction getE() {
        return this.e;
    }

    public NDDOSolution getS() {
        return s;
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

    public double[][] getTotalGradients() {
        return totalGradients;
    }
}
