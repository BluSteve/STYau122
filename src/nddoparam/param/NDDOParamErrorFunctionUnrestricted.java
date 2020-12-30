package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDODerivative;
import nddoparam.NDDOSolution;

public class NDDOParamErrorFunctionUnrestricted extends NDDOParamErrorFunction {
    public NDDOParamErrorFunctionUnrestricted(NDDOAtom[] atoms, NDDOSolution soln, double refHeat) {
        super(atoms, soln, refHeat);
    }

    @Override
    protected double getGradient(int i, int j) {
        return NDDODerivative.gradientUnrestricted(expAtoms, expSoln.alphaDensity(), expSoln.betaDensity(), i, j);
    }

    @Override
    protected double getDeriv(double[] coeff, int atom2) {
        double deriv = coeff[0] * NDDODerivative.gradientUnrestricted(expAtoms, expSoln.alphaDensity(), expSoln.betaDensity(), atom2, 0)
                + coeff[1] * NDDODerivative.gradientUnrestricted(expAtoms, expSoln.alphaDensity(), expSoln.betaDensity(), atom2, 1)
                + coeff[2] * NDDODerivative.gradientUnrestricted(expAtoms, expSoln.alphaDensity(), expSoln.betaDensity(), atom2, 2);
        return 1E-13 * Math.round(deriv * 1E13);
    }
}
