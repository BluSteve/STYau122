package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.NDDOSolution;
import scf.Utils;

public abstract class ParamGradient {
    protected NDDOAtom[] atoms, perturbed;
    protected int Z, paramNum;
    protected NDDOSolution s, sprime;
    protected ParamErrorFunction e, eprime;

    public ParamGradient(NDDOAtom[] atoms, int Z, int paramNum) {
        this.Z = Z;
        this.paramNum = paramNum;
        this.atoms = atoms;
        perturbed = Utils.perturbAtomParams(atoms, paramNum, Z);
    }

    public void addDipoleError(double ref) {
        if (eprime != null) {
            eprime.addDipoleError(ref);
        }
        e.addDipoleError(ref);
    }

    public void addIEError(double ref) {
        if (eprime != null) {
            eprime.addIEError(ref);
        }
        e.addIEError(ref);
    }

    public void addGeomError() {
        if (eprime != null) {
            eprime.addGeomError();
        }
        e.addGeomError();
    }

    public void addBondError(int atom1, int atom2, double ref) {
        if (eprime != null) {
            eprime.addBondError(atom1, atom2, ref);
        }
        e.addBondError(atom1, atom2, ref);
    }

    public void addAngleError(int atom1, int atom2, int atom3, double ref) {
        if (eprime != null) {
            eprime.addAngleError(atom1, atom2, atom3, ref);
        }
        e.addAngleError(atom1, atom2, atom3, ref);
    }

    public double gradient() {
        if (eprime != null) {
            return (eprime.getTotalError() - e.getTotalError()) / Utils.lambda;
        } else return 0;
    }

    public abstract void constructErrors(double refHeat);

    public abstract void createExpGeom(NDDOAtom[] expAtoms, NDDOSolution expSoln);

    public double getEDipole() {
        return s.dipole;
    }

    public double getEPrimeDipole() {
        if (eprime != null) {
            return eprime.soln.dipole;
        }
        return s.dipole;
    }

    public double getEHf() {
        return s.hf;
    }

    public double getEPrimeHf() {
        if (eprime != null) {
            return eprime.soln.hf;
        }
        return s.hf;
    }

    public double getEGradient() {
        return e.geomGradient;
    }

    public double getEPrimeGradient() {
        if (eprime != null) {
            return eprime.geomGradient;
        }
        return e.geomGradient;
    }

    public double getEHomo() {
        return s.homo;
    }

    public double getEPrimeHomo() {
        if (eprime != null) {
            return eprime.soln.homo;
        }
        return s.homo;
    }

    public ParamErrorFunction getE() {
        return e;
    }

    public ParamErrorFunction getEPrime() {
        return eprime;
    }
}


