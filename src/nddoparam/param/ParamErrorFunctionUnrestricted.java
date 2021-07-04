package nddoparam.param;

import nddoparam.NDDODerivative;
import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionUnrestricted;

public class ParamErrorFunctionUnrestricted extends ParamErrorFunction {
    public ParamErrorFunctionUnrestricted(NDDOSolution soln, double refHeat) {
        super(soln, refHeat);
    }

    @Override
    protected double getGradient(int i, int j) {
        return NDDODerivative.grad(expAtoms, (NDDOSolutionUnrestricted) expSoln, i, j);
    }

    @Override
    protected double getDeriv(double[] coeff, int atom2) {
        double deriv = coeff[0] * NDDODerivative.grad(expAtoms, (NDDOSolutionUnrestricted) expSoln, atom2, 0)
                + coeff[1] * NDDODerivative.grad(expAtoms, (NDDOSolutionUnrestricted) expSoln, atom2, 1)
                + coeff[2] * NDDODerivative.grad(expAtoms, (NDDOSolutionUnrestricted) expSoln, atom2, 2);
        return 1E-13 * Math.round(deriv * 1E13);
    }
}
