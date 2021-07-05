package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimizationRestricted;
import nddoparam.NDDOSolutionRestricted;
import nddoparam.param.ParamGradientRestricted;
import nddoparam.param.ParamHessianRestricted;

public class MoleculeRunRestricted extends MoleculeRun {
    public MoleculeRunRestricted(NDDOAtom[] atoms2, int charge, NDDOAtom[] expGeom, double[] datum, boolean runHessian, String kind, int[] atomTypes) {
        super(atoms2, charge, expGeom, datum, runHessian, kind, atomTypes, 1);
        opt = new NDDOGeometryOptimizationRestricted(atoms, charge); // geometry is optimized upon initialization

        super.newGeomCoords = "RHF\n" + "CHARGE=" + charge + "\nMULT=1\n";
        super.generateGeomCoords();
        if (expGeom != null) expSolution = new NDDOSolutionRestricted(this.expGeom, charge);

        if (runHessian) {
            runHessian();
        } else {
            hessianStr = "";
        }
        runGradient();

        outputErrorFunction();
    }

    protected void getH(int Z1, int paramNum1, int Z2, int paramNum2) {
        h = new ParamHessianRestricted(atoms, charge, Z1, paramNum1, Z2, paramNum2, (NDDOSolutionRestricted) opt.s);
    }

}
