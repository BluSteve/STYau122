package mndoparam.param;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDOParams;
import mndoparam.mndo.MNDOSolutionUnrestricted;

public class MNDOParamHessianUnrestricted {

    private static final double lambda = 1E-7;

    public MNDOParamGradientUnrestricted g, gprime;

    private MNDOAtom[] atoms, perturbed;

    public String str;

    private int Z1;

    private int paramnum1;

    public MNDOParamHessianUnrestricted(MNDOAtom[] atoms, int charge, int mult, int Z1, int paramnum1, int Z2, int paramnum2, MNDOSolutionUnrestricted s) {

        this.Z1 = Z1;

        this.paramnum1 = paramnum1;

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

        MNDOSolutionUnrestricted sprime = new MNDOSolutionUnrestricted(perturbed, charge, mult);

        g = new MNDOParamGradientUnrestricted(atoms, charge, mult, Z2, paramnum2, s);

        gprime = new MNDOParamGradientUnrestricted(perturbed, charge, mult, Z2, paramnum2, sprime);

        this.str = "";

        switch (Z1) {
            case 1:
                str += "H ";
                break;
            case 6:
                str += "C ";
                break;
            case 7:
                str += "N ";
                break;
        }

        switch (paramnum1) {
            case 0:
                str += "ALPHA, ";
                break;
            case 1:
                str += "BETAS, ";
                break;
            case 2:
                str += "BETAP, ";
                break;
            case 3:
                str += "USS, ";
                break;
            case 4:
                str += "UPP, ";
                break;
            case 5:
                str += "ZETAS, ";
                break;
            case 6:
                str += "ZETAP, ";
                break;
            case 7:
                str += "EISOL, ";
                break;
            default:
                System.exit(0);

        }

        switch (Z2) {
            case 1:
                str += "H ";
                break;
            case 6:
                str += "C ";
                break;
            case 7:
                str += "N ";
                break;
        }

        switch (paramnum2) {
            case 0:
                str += "ALPHA: ";
                break;
            case 1:
                str += "BETAS: ";
                break;
            case 2:
                str += "BETAP: ";
                break;
            case 3:
                str += "USS: ";
                break;
            case 4:
                str += "UPP: ";
                break;
            case 5:
                str += "ZETAS: ";
                break;
            case 6:
                str += "ZETAP: ";
                break;
            case 7:
                str += "EISOL: ";
                break;
            default:
                System.exit(0);

        }

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

    public void createExpGeom(MNDOAtom[] expatoms, MNDOSolutionUnrestricted expsoln) {
        g.createExpGeom(expatoms, expsoln);
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
        gprime.createExpGeom(perturbed, new MNDOSolutionUnrestricted(perturbed, expsoln.charge, expsoln.multiplicity));

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