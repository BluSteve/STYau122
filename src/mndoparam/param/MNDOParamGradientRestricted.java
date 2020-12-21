package mndoparam.param;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDOSolution;
import mndoparam.mndo.MNDOSolutionRestricted;
import scf.Utils;

public class MNDOParamGradientRestricted extends MNDOParamGradient {
    public MNDOParamGradientRestricted(MNDOAtom[] atoms, int charge, int Z, int paramNum, MNDOSolutionRestricted s) {
        super(atoms, Z, paramNum);
        this.s = s;
        sprime = new MNDOSolutionRestricted(perturbed, charge);
    }

    public void constructErrors(double refHeat) {
        e = new MNDOParamErrorFunctionRestricted(atoms, s, refHeat);
        eprime = new MNDOParamErrorFunctionRestricted(perturbed, sprime, refHeat);
    }

    @Override
    public void createExpGeom(MNDOAtom[] expAtoms, MNDOSolution expSoln) {
        e.createExpGeom(expAtoms, expSoln);
        MNDOAtom[] perturbedExp = Utils.perturbAtoms(expAtoms, paramNum, Z);
        eprime.createExpGeom(perturbedExp, new MNDOSolutionRestricted(perturbedExp, expSoln.charge));
    }
}
