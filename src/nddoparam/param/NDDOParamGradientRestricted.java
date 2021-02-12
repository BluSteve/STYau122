package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;
import scf.Utils;

public class NDDOParamGradientRestricted extends NDDOParamGradient {
    public NDDOParamGradientRestricted(NDDOAtom[] atoms, int charge, int Z, int paramNum, NDDOSolutionRestricted s) {
        super(atoms, Z, paramNum);
        this.s = s;
        sprime = new NDDOSolutionRestricted(perturbed, charge);
    }

    public void constructErrors(double refHeat) {
        if (Utils.containsZ(atoms, Z)) {
            e = new NDDOParamErrorFunctionRestricted(atoms, s, refHeat);
            eprime = new NDDOParamErrorFunctionRestricted(perturbed, sprime, refHeat);
        } else {
            e = new NDDOParamErrorFunctionRestricted(atoms, s, refHeat);
        }
    }

    @Override
    public void createExpGeom(NDDOAtom[] expAtoms, NDDOSolution expSoln) {
        if (eprime != null) {
            NDDOAtom[] perturbedExp = Utils.perturbAtomParams(expAtoms, paramNum, Z);
            eprime.createExpGeom(perturbedExp, new NDDOSolutionRestricted(perturbedExp, expSoln.charge));
        }
        e.createExpGeom(expAtoms, expSoln);
    }
}
