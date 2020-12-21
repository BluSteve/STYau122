package mndoparam.param;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDOSolution;
import mndoparam.mndo.MNDOSolutionRestricted;
import scf.Utils;

public class MNDOParamHessianRestricted extends MNDOParamHessian {
    public MNDOParamHessianRestricted(MNDOAtom[] atoms, int charge, int Z1, int paramNum1, int Z2, int paramNum2, MNDOSolutionRestricted s) {
        super(atoms, Z1, paramNum1);
        MNDOSolutionRestricted sPrime = new MNDOSolutionRestricted(perturbed, charge);
        g = new MNDOParamGradientRestricted(atoms, charge, Z2, paramNum2, s);
        gprime = new MNDOParamGradientRestricted(perturbed, charge, Z2, paramNum2, sPrime);
        //System.err.println("initialization complete");
    }

    public void createExpGeom(MNDOAtom[] expAtoms, MNDOSolution expSolution) {
        //System.err.println("creating expgeom");
        g.createExpGeom(expAtoms, expSolution);
        MNDOAtom[] perturbed = Utils.perturbAtoms(expAtoms, paramNum1, Z1);
        gprime.createExpGeom(perturbed, new MNDOSolutionRestricted(perturbed, expSolution.charge));
        //System.err.println("creation complete");
    }
}
