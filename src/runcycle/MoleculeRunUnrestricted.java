package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimizationUnrestricted;
import nddoparam.NDDOSolutionUnrestricted;
import nddoparam.param.ParamGradientUnrestricted;
import nddoparam.param.ParamHessianUnrestricted;

public class MoleculeRunUnrestricted extends MoleculeRun {
    public MoleculeRunUnrestricted(NDDOAtom[] atoms2, int charge, int mult, NDDOAtom[] expGeom, double[] datum, boolean runHessian, String kind, int[] atomTypes) {
        super(atoms2, charge, expGeom, datum, runHessian, kind, atomTypes, mult);
        opt = new NDDOGeometryOptimizationUnrestricted(atoms, charge, mult);

        newGeomCoords = "UHF\n" + "CHARGE=" + charge + "\nMULT=" + mult + "\n";
        generateGeomCoords();
        if (this.expGeom != null) expSolution = new NDDOSolutionUnrestricted(this.expGeom, charge, mult);

        if (runHessian) {
            runHessian();
        } else {
            hessianStr = "";
        }
        runGradient();

        outputErrorFunction();
    }

    protected void getG(int Z, int numit, int mult) {
        // TODO uncomment this
//        g = new ParamGradientUnrestricted(atoms, charge, mult, Z, numit, (NDDOSolutionUnrestricted) opt.s);

    }

    protected void getH(int Z1, int param1, int Z2, int param2) {
        h = new ParamHessianUnrestricted(atoms, charge, mult, Z1, param1, Z2, param2, (NDDOSolutionUnrestricted) opt.s);
    }

}