package mndoparam.param;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDOSolution;
import mndoparam.mndo.MNDOSolutionUnrestricted;
import scf.Utils;

public class MNDOParamHessianUnrestricted extends MNDOParamHessian {
    public MNDOParamHessianUnrestricted(MNDOAtom[] atoms, int charge, int mult, int Z1, int paramNum1, int Z2, int paramNum2, MNDOSolutionUnrestricted s) {
        super(atoms, Z1, paramNum1);
        MNDOSolutionUnrestricted sPrime = new MNDOSolutionUnrestricted(perturbed, charge, mult);
        g = new MNDOParamGradientUnrestricted(atoms, charge, mult, Z2, paramNum2, s);
        gprime = new MNDOParamGradientUnrestricted(perturbed, charge, mult, Z2, paramNum2, sPrime);
        //System.err.println("initialization complete");
    }

    public void createExpGeom(MNDOAtom[] expAtoms, MNDOSolution expSolution) {
        g.createExpGeom(expAtoms, expSolution);
        //System.err.println("creating expgeom");
        MNDOAtom[] perturbed = Utils.perturbAtoms(expAtoms, paramNum1, Z1);
        gprime.createExpGeom(perturbed, new MNDOSolutionUnrestricted(perturbed, expSolution.charge, expSolution.multiplicity));
        //System.err.println("creation complete");
    }


}