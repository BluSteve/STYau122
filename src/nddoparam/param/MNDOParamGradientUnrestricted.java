package nddoparam.param;

import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOSolution;
import nddoparam.mndo.MNDOSolutionUnrestricted;
import scf.Utils;

public class MNDOParamGradientUnrestricted extends MNDOParamGradient {
    public MNDOParamGradientUnrestricted(MNDOAtom[] atoms, int charge, int mult, int Z, int paramNum, MNDOSolutionUnrestricted s) {
        super(atoms, Z, paramNum);
        this.s = s;
        sprime = new MNDOSolutionUnrestricted(perturbed, charge, mult);
    }

    public void constructErrors(double refHeat) {
        e = new MNDOParamErrorFunctionUnrestricted(atoms, s, refHeat);
        eprime = new MNDOParamErrorFunctionUnrestricted(perturbed, sprime, refHeat);
    }

    // TODO AbstractMNDOSolution will probably not work
    @Override
    public void createExpGeom(MNDOAtom[] expAtoms, MNDOSolution expSoln) {
        e.createExpGeom(expAtoms, expSoln);
        MNDOAtom[] perturbedExp = Utils.perturbAtoms(expAtoms, paramNum, Z);
        eprime.createExpGeom(perturbedExp, new MNDOSolutionUnrestricted(perturbedExp, expSoln.charge, expSoln.multiplicity));
    }
}
