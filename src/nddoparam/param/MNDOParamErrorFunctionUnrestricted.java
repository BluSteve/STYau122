package nddoparam.param;

import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDODerivative;
import nddoparam.mndo.MNDOSolution;

public class MNDOParamErrorFunctionUnrestricted extends MNDOParamErrorFunction {
    public MNDOParamErrorFunctionUnrestricted(MNDOAtom[] atoms, MNDOSolution soln, double refHeat) {
        super(atoms, soln, refHeat);
    }

    @Override
    protected double getGradient(int i, int j) {
        return MNDODerivative.gradientUnrestricted(expAtoms, expSoln.alphaDensity(), expSoln.betaDensity(), i, j);
    }

    @Override
    protected double getDeriv(double[] coeff, int atom2) {
        double deriv = coeff[0] * MNDODerivative.gradientUnrestricted(expAtoms, expSoln.alphaDensity(), expSoln.betaDensity(), atom2, 0)
                + coeff[1] * MNDODerivative.gradientUnrestricted(expAtoms, expSoln.alphaDensity(), expSoln.betaDensity(), atom2, 1)
                + coeff[2] * MNDODerivative.gradientUnrestricted(expAtoms, expSoln.alphaDensity(), expSoln.betaDensity(), atom2, 2);
        return 1E-13 * Math.round(deriv * 1E13);
    }
}
