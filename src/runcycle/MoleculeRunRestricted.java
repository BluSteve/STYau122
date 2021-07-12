package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimizationRestricted;
import nddoparam.NDDOSolutionRestricted;
import nddoparam.param.ParamGradientRestricted;
import nddoparam.param.ParamHessianRestricted;

public class MoleculeRunRestricted extends MoleculeRun {
    public MoleculeRunRestricted(NDDOAtom[] atoms2, int charge, NDDOAtom[] expGeom, double[] datum, boolean isRunHessian, String kind, int[] atomTypes) {
        super(atoms2, charge, expGeom, datum, isRunHessian, kind, atomTypes, 1);

        opt = new NDDOGeometryOptimizationRestricted(atoms, charge); // ~ 74 ms for CH4 type d with expsoln

        newGeomCoords = "RHF\n" + "CHARGE=" + charge + "\nMULT=1\n";
        generateGeomCoords();
        if (expGeom != null) expSolution = new NDDOSolutionRestricted(this.expGeom, charge); // ~ 60 ms

        routine(); // ~700-800 ms
    }

    @Override
    protected void constructG() {
        g = new ParamGradientRestricted((NDDOSolutionRestricted) opt.s, kind, datum, (NDDOSolutionRestricted) expSolution, true);
    }

    @Override
    protected void constructH() {
        h = new ParamHessianRestricted((ParamGradientRestricted) g, true);
    }
}
