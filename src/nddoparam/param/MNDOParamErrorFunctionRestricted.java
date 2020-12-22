package nddoparam.param;

import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDODerivative;
import nddoparam.mndo.MNDOSolution;

public class MNDOParamErrorFunctionRestricted extends MNDOParamErrorFunction {
    public MNDOParamErrorFunctionRestricted(MNDOAtom[] atoms, MNDOSolution soln, double refHeat) {
        super(atoms, soln, refHeat);
    }


    @Override
    protected double getGradient(int i, int j) {
        return MNDODerivative.gradient(expAtoms, expSoln.densityMatrix(), i, j);
    }

    @Override
    protected double getDeriv(double[] coeff, int atom2) {
        double deriv = coeff[0] * MNDODerivative.gradient(expAtoms, expSoln.densityMatrix(), atom2, 0)
                + coeff[1] * MNDODerivative.gradient(expAtoms, expSoln.densityMatrix(), atom2, 1)
                + coeff[2] * MNDODerivative.gradient(expAtoms, expSoln.densityMatrix(), atom2, 2);
        return 1E-13 * Math.round(deriv * 1E13);
    }

}
