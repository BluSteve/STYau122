package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;
import scf.Utils;

public class NDDOParamHessianRestricted extends NDDOParamHessian {
    public NDDOParamHessianRestricted(NDDOAtom[] atoms, int charge, int Z1, int paramNum1, int Z2, int paramNum2, NDDOSolutionRestricted s) {
        super(atoms, Z1, paramNum1);
        NDDOSolutionRestricted sPrime = new NDDOSolutionRestricted(perturbed, charge);
        g = new NDDOParamGradientRestricted(atoms, charge, Z2, paramNum2, s);
        gprime = new NDDOParamGradientRestricted(perturbed, charge, Z2, paramNum2, sPrime);
        //System.err.println("initialization complete");
    }

    public void createExpGeom(NDDOAtom[] expAtoms, NDDOSolution expSolution) {
        //System.err.println("creating expgeom");
        g.createExpGeom(expAtoms, expSolution);
        NDDOAtom[] perturbed = Utils.perturbAtoms(expAtoms, paramNum1, Z1);
        gprime.createExpGeom(perturbed, new NDDOSolutionRestricted(perturbed, expSolution.charge));
        //System.err.println("creation complete");
    }
}
