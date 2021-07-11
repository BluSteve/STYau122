package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDOSolution;
import scf.Utils;

public abstract class ParamHessian implements ErrorGettable{
    public ParamGradient g, gprime;
    protected NDDOAtom[] atoms, perturbed;
    protected int paramNum1, Z1;

    public ParamHessian(NDDOAtom[] atoms, int Z1, int paramNum1) {
        this.Z1 = Z1;
        this.paramNum1 = paramNum1;
        this.atoms = atoms;
        perturbed = Utils.perturbAtomParams(atoms, Z1, paramNum1);
    }

    public void constructErrors(double refHeat) {
        g.constructErrors(refHeat);
        gprime.constructErrors(refHeat);
    }

    public void addDipoleError(double ref) {
        g.addDipoleError(ref);
        gprime.addDipoleError(ref);
    }

    public void addIEError(double ref) {
        g.addIEError(ref);
        gprime.addIEError(ref);
    }

    public void addGeomError() {
        g.addGeomError();
        gprime.addGeomError();
    }

    public void addBondError(int atom1, int atom2, double ref) {
        g.addBondError(atom1, atom2, ref);
        gprime.addBondError(atom1, atom2, ref);
    }

    public void addAngleError(int atom1, int atom2, int atom3, double ref) {
        g.addAngleError(atom1, atom2, atom3, ref);
        gprime.addAngleError(atom1, atom2, atom3, ref);
    }

    public abstract void createExpGeom(NDDOAtom[] expAtoms, NDDOSolution expSoln);

    public double hessian() {
        return (gprime.gradient() - g.gradient()) / Utils.lambda;
    }

}
