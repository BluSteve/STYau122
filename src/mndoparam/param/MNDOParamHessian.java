package mndoparam.param;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDOParams;
import mndoparam.mndo.MNDOSolution;

public class MNDOParamHessian {

    private static final double lambda = 1E-7;

    public MNDOParamGradient g, gprime;

    private MNDOAtom[] atoms, perturbed;

    public String str;

    private int paramnum1, Z1;

    public MNDOParamHessian(MNDOAtom[] atoms, int charge, int Z1, int paramnum1, int Z2, int paramnum2, MNDOSolution s) {

        System.err.println("initializing Hessian");

        this.paramnum1 = paramnum1;

        this.Z1 = Z1;

        this.atoms = atoms;

        perturbed = new MNDOAtom[atoms.length];

        for (int i = 0; i < atoms.length; i++) {

            perturbed[i] = new MNDOAtom(atoms[i]);

            if (atoms[i].getAtomProperties().getZ() == Z1) {
                MNDOParams params = atoms[i].getParams();

                params.modifyParam(paramnum1, lambda);

                perturbed[i] = new MNDOAtom(atoms[i], params);
            }
        }

        MNDOSolution sprime = new MNDOSolution(perturbed, charge);

        g = new MNDOParamGradient(atoms, charge, Z2, paramnum2, s);

        gprime = new MNDOParamGradient(perturbed, charge, Z2, paramnum2, sprime);

        System.err.println("initialization complete");
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

    public void createExpGeom(MNDOAtom[] expatoms, MNDOSolution expsoln) {

        System.err.println("creating expgeom");

        g.createExpGeom(expatoms, expsoln);

        MNDOAtom[] perturbed = new MNDOAtom[expatoms.length];

        for (int i = 0; i < expatoms.length; i++) {

            perturbed[i] = new MNDOAtom(expatoms[i]);

            if (expatoms[i].getAtomProperties().getZ() == Z1) {
                MNDOParams params = expatoms[i].getParams();
                params.modifyParam(paramnum1, lambda);
                perturbed[i] = new MNDOAtom(expatoms[i], params);
            }
        }
        gprime.createExpGeom(perturbed, new MNDOSolution(perturbed, expsoln.charge));

        System.err.println("creation complete");
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

    public double hessian() {
        return (gprime.gradient() - g.gradient()) / lambda;
    }

}
