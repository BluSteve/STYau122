package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionUnrestricted;
import scf.Utils;

public class ParamGradientUnrestricted extends ParamGradient {
    public ParamGradientUnrestricted(NDDOAtom[] atoms, int charge, int mult, int Z, int paramNum, NDDOSolutionUnrestricted s) {
        super(atoms, Z, paramNum);
        this.s = s;
        sprime = new NDDOSolutionUnrestricted(perturbed, charge, mult);
    }

    public void constructErrors(double refHeat) {
        if (Utils.containsZ(atoms, Z)) {
            e = new ParamErrorFunctionUnrestricted(atoms, s, refHeat);
            eprime = new ParamErrorFunctionUnrestricted(perturbed, sprime, refHeat);
        }
        else {
            e = new ParamErrorFunctionUnrestricted(atoms, s, refHeat);
        }
    }

    @Override
    public void createExpGeom(NDDOAtom[] expAtoms, NDDOSolution expSoln) {
        if (eprime != null) {
            NDDOAtom[] perturbedExp = Utils.perturbAtomParams(expAtoms, paramNum, Z);
            eprime.createExpGeom(perturbedExp, new NDDOSolutionUnrestricted(perturbedExp, expSoln.charge, expSoln.multiplicity));
        }
        e.createExpGeom(expAtoms, expSoln);
    }
}
