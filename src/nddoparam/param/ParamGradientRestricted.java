package nddoparam.param;

import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;

public class ParamGradientRestricted extends ParamGradientAnalytical {
    public ParamGradientRestricted(NDDOSolution s, String kind) {
        this(s, kind, 0);
    }

    public ParamGradientRestricted(NDDOSolution s, String kind, int charge) {
        super(s, kind);
        this.charge = charge;
    }

    @Override
    public void constructErrors(double refHeat) {
        this.e = new ParamErrorFunctionRestricted(s, refHeat);
    }

    // gets only HFDerivs
    @Override
    protected void computeHFDerivs() {
        for (int Z : s.getUniqueZs()) {
            for (int paramNum : s.getNeededParams()[Z]) {
                this.HFDerivs[Z][paramNum] = ParamDerivative.HfDeriv((NDDOSolutionRestricted) this.s, Z, paramNum);
            }
        }
    }

    // gets HFDerivs and dipoleDerivs
    @Override
    protected void computeLimitedDerivs() {
        for (int Z : s.getUniqueZs()) {
            computeBatchedDerivs(Z);

            for (int paramNum : s.getNeededParams()[Z]) {
                computeDipoleDeriv(Z, paramNum);
            }
        }
    }

    // gets HFDerivs, dipoleDerivs and IEDerivs;
    @Override
    protected void computeComplementaryDerivs() {
        for (int Z : s.getUniqueZs()) {
            computeBatchedDerivs(Z);

            for (int paramNum : s.getNeededParams()[Z]) {
                computeDipoleDeriv(Z, paramNum);
                computeIEDeriv(Z, paramNum);
            }
        }
    }

    @Override
    protected void computeBatchedDerivs(int Z) {
        this.staticDerivs[Z] = ParamDerivative.MNDOStaticMatrixDeriv((NDDOSolutionRestricted) this.s, Z);
        this.xLimited[Z] = ParamDerivative.xarraylimited((NDDOSolutionRestricted) this.s,
                this.staticDerivs[Z][1]);
    }

    @Override
    protected void computeDipoleDeriv(int Z, int paramNum) {
        this.HFDerivs[Z][paramNum] = ParamDerivative.MNDOHfDeriv((NDDOSolutionRestricted) this.s,
                this.staticDerivs[Z][0][paramNum], this.staticDerivs[Z][1][paramNum]);
        this.densityDerivs[Z][paramNum] = ParamDerivative.densityDerivativeLimited((NDDOSolutionRestricted) this.s,
                this.xLimited[Z][paramNum]);
        this.dipoleDerivs[Z][paramNum] = ParamDerivative.MNDODipoleDeriv((NDDOSolutionRestricted) this.s, densityDerivs[Z][paramNum], Z,
                paramNum);
    }

    @Override
    protected void computeIEDeriv(int Z, int paramNum) {
        this.responseDerivs[Z][paramNum] = ParamDerivative.responseMatrix((NDDOSolutionRestricted)
                this.s, this.densityDerivs[Z][paramNum]);
        this.fockDerivs[Z][paramNum] = this.staticDerivs[Z][1][paramNum].add(this.responseDerivs[Z][paramNum]);
        this.xComplementary[Z][paramNum] = ParamDerivative.xArrayComplementary((NDDOSolutionRestricted) this.s,
                this.fockDerivs[Z][paramNum]);
        this.xForIE[Z][paramNum] = ParamDerivative.xarrayForIE((NDDOSolutionRestricted) this.s,
                xLimited[Z][paramNum], xComplementary[Z][paramNum]);
        this.coeffDerivs[Z][paramNum] = ParamDerivative.HOMOcoefficientDerivativeComplementary(this.xForIE[Z][paramNum],
                (NDDOSolutionRestricted) this.s);
        this.IEDerivs[Z][paramNum] = ParamDerivative.MNDOIEDeriv((NDDOSolutionRestricted) this.s,
                this.coeffDerivs[Z][paramNum], this.fockDerivs[Z][paramNum]);
    }
}
