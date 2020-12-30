package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDODerivative;
import nddoparam.NDDOSolution;

public class NDDOParamErrorFunctionRestricted extends NDDOParamErrorFunction {
    public NDDOParamErrorFunctionRestricted(NDDOAtom[] atoms, NDDOSolution soln, double refHeat) {
        super(atoms, soln, refHeat);
    }

    @Override
    protected double getGradient(int i, int j) {
        return NDDODerivative.gradient(expAtoms, expSoln.densityMatrix(), i, j);
    }

    @Override
    protected double getDeriv(double[] coeff, int atom2) {
        double deriv = coeff[0] * NDDODerivative.gradient(expAtoms, expSoln.densityMatrix(), atom2, 0)
                + coeff[1] * NDDODerivative.gradient(expAtoms, expSoln.densityMatrix(), atom2, 1)
                + coeff[2] * NDDODerivative.gradient(expAtoms, expSoln.densityMatrix(), atom2, 2);
        return 1E-13 * Math.round(deriv * 1E13);
    }

}
