package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDODerivative;
import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;

public class NDDOParamErrorFunctionRestricted extends NDDOParamErrorFunction {
    public NDDOParamErrorFunctionRestricted(NDDOAtom[] atoms, NDDOSolution soln, double refHeat) {
        super(atoms, soln, refHeat);
    }

    @Override
    protected double getGradient(int i, int j) {
        return NDDODerivative.gradient(expAtoms, (NDDOSolutionRestricted) expSoln, i, j);
    }

    @Override
    protected double getDeriv(double[] coeff, int atom2) {
        double deriv = coeff[0] * NDDODerivative.gradient(expAtoms, (NDDOSolutionRestricted) expSoln, atom2, 0)
                + coeff[1] * NDDODerivative.gradient(expAtoms, (NDDOSolutionRestricted) expSoln, atom2, 1)
                + coeff[2] * NDDODerivative.gradient(expAtoms, (NDDOSolutionRestricted) expSoln, atom2, 2);
        return 1E-13 * Math.round(deriv * 1E13);
    }

}
