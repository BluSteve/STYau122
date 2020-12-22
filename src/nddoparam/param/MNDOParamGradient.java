package nddoparam.param;

import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOSolution;
import scf.Utils;

public abstract class MNDOParamGradient {
    protected MNDOAtom[] atoms, perturbed;
    protected int Z, paramNum;
    public MNDOSolution s, sprime;

    public MNDOParamErrorFunction e, eprime;

    // TODO what is charge and Z?
    public MNDOParamGradient(MNDOAtom[] atoms, int Z, int paramNum) {
        this.Z = Z;
        this.paramNum = paramNum;
        this.atoms = atoms;
        perturbed = Utils.perturbAtoms(atoms, paramNum, Z);
    }

    public void addDipoleError(double ref) {
        e.AddDipoleError(ref);
        eprime.AddDipoleError(ref);
    }

    public void addIEError(double ref) {
        e.AddIEError(ref);
        eprime.AddIEError(ref);
    }

    public void addGeomError() {
        e.addGeomError();
        eprime.addGeomError();
    }

    public void addBondError(int atom1, int atom2, double ref) {
        e.addBondError(atom1, atom2, ref);
        eprime.addBondError(atom1, atom2, ref);
    }

    public void addAngleError(int atom1, int atom2, int atom3, double ref) {
        e.addAngleError(atom1, atom2, atom3, ref);
        eprime.addAngleError(atom1, atom2, atom3, ref);
    }

    public double gradient() {
        return (eprime.constructErrorFunction() - e.constructErrorFunction()) / Utils.lambda;
    }

    public abstract void constructErrors(double refHeat);

    public abstract void createExpGeom(MNDOAtom[] expAtoms, MNDOSolution expSoln);
}
