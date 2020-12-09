package mndoparam.param;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDOGeometryOptimization;
import mndoparam.mndo.MNDOParams;
import mndoparam.mndo.MNDOSolution;

public class MNDOParamGradient {

    public static final double lambda = 1E-7;

    public MNDOSolution s;

    private MNDOSolution sprime;

    public MNDOParamErrorFunction e, eprime;

    private MNDOAtom[] atoms;

    public MNDOAtom[] perturbed;

    private int Z, paramNum;

    // TODO reduce number of constructors
    public MNDOParamGradient(MNDOAtom[] atoms, int charge, int Z, int paramNum, MNDOGeometryOptimization opt) {
        this.Z = Z;
        this.paramNum = paramNum;

        this.atoms = atoms;

        this.s = opt.s;

        perturbed = new MNDOAtom[atoms.length];

        for (int i = 0; i < atoms.length; i++) {

            perturbed[i] = new MNDOAtom(atoms[i]);

            if (atoms[i].getAtomProperties().getZ() == Z) {
                MNDOParams params = atoms[i].getParams();

                params.modifyParam(paramNum, lambda);

                perturbed[i] = new MNDOAtom(atoms[i], params);
            }
        }

        sprime = new MNDOSolution(perturbed, charge);
    }

    public MNDOParamGradient(MNDOAtom[] atoms, int charge, int Z, int paramNum, MNDOSolution s) {
        this.Z = Z;
        this.paramNum = paramNum;

        this.atoms = atoms;

        this.s = s;

        perturbed = new MNDOAtom[atoms.length];

        for (int i = 0; i < atoms.length; i++) {

            perturbed[i] = new MNDOAtom(atoms[i]);

            if (atoms[i].getAtomProperties().getZ() == Z) {
                MNDOParams params = atoms[i].getParams();

                params.modifyParam(paramNum, lambda);

                perturbed[i] = new MNDOAtom(atoms[i], params);
            }
        }

        sprime = new MNDOSolution(perturbed, charge);
    }

    public void constructErrors(double refHeat) {
        e = new MNDOParamErrorFunction(atoms, s, refHeat);
        eprime = new MNDOParamErrorFunction(perturbed, sprime, refHeat);
    }

    public void addDipoleError(double ref) {

        e.AddDipoleError(ref);
        eprime.AddDipoleError(ref);
    }

    public void addIEError(double ref) {
        e.AddIEError(ref);
        eprime.AddIEError(ref);
    }

    public void createExpGeom(MNDOAtom[] expatoms, MNDOSolution expsoln) {
        e.createExpGeom(expatoms, expsoln);

        MNDOAtom[] perturbed = new MNDOAtom[expatoms.length];

        for (int i = 0; i < expatoms.length; i++) {

            perturbed[i] = new MNDOAtom(expatoms[i]);

            //System.err.println (Z);

            if (expatoms[i].getAtomProperties().getZ() == Z) {
                MNDOParams params = expatoms[i].getParams();
                params.modifyParam(paramNum, lambda);
                perturbed[i] = new MNDOAtom(expatoms[i], params);
            }
        }
        eprime.createExpGeom(perturbed, new MNDOSolution(perturbed, expsoln.charge));
    }

    public void addGeomError() {
        e.AddGeomError();
        eprime.AddGeomError();
    }

    public void addBondError(int atom1, int atom2, double ref) {
        e.AddBondError(atom1, atom2, ref);
        eprime.AddBondError(atom1, atom2, ref);
    }

    public void addAngleError(int atom1, int atom2, int atom3, double ref) {
        e.AddAngleError(atom1, atom2, atom3, ref);
        eprime.AddAngleError(atom1, atom2, atom3, ref);
    }

    public double gradient() {
        return (eprime.constructErrorFunction() - e.constructErrorFunction()) / lambda;
    }

}
