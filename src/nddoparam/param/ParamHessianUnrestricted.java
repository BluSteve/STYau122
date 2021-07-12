package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionUnrestricted;
import scf.Utils;

public class ParamHessianUnrestricted extends ParamHessian {
    public ParamHessianUnrestricted(NDDOAtom[] atoms, int charge, int mult, int Z1, int paramNum1, int Z2, int paramNum2, NDDOSolutionUnrestricted s) {
        super(atoms, Z1, paramNum1);
        NDDOSolutionUnrestricted sPrime = new NDDOSolutionUnrestricted(perturbed, charge, mult);
        g = new ParamGradientUnrestricted(atoms, charge, mult, Z2, paramNum2, s);
        gprime = new ParamGradientUnrestricted(perturbed, charge, mult, Z2, paramNum2, sPrime);
    }

    public void createExpGeom(NDDOAtom[] expAtoms, NDDOSolution expSolution) {
        g.createExpGeom(expAtoms, expSolution);
        NDDOAtom[] perturbed = Utils.perturbAtomParams(expAtoms, Z1, paramNum1);
        gprime.createExpGeom(perturbed, new NDDOSolutionUnrestricted(perturbed, expSolution.charge, expSolution.multiplicity));
    }


    @Override
    public ParamErrorFunction getE() {
        return null;
    }
}