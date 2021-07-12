package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimizationUnrestricted;
import nddoparam.NDDOSolutionUnrestricted;
import nddoparam.param.*;

public class MoleculeRunUnrestricted extends MoleculeRun {
    public MoleculeRunUnrestricted(NDDOAtom[] atoms2, int charge, int mult, NDDOAtom[] expGeom, double[] datum, boolean isRunHessian, String kind, int[] atomTypes) {
        super(atoms2, charge, expGeom, datum, isRunHessian, kind, atomTypes, mult);
        opt = new NDDOGeometryOptimizationUnrestricted(atoms, charge, mult);

        newGeomCoords = "UHF\n" + "CHARGE=" + charge + "\nMULT=" + mult + "\n";
        generateGeomCoords();
        if (this.expGeom != null) expSolution = new NDDOSolutionUnrestricted(this.expGeom, charge, mult);

        if (isRunHessian) {
            runHessian();
        } else {
            hessianStr = "";
        }
        runGradient();

        outputErrorFunction();
    }

    @Override
    protected void constructG() {
        g = new ParamGradientUnrestricted2((NDDOSolutionUnrestricted) opt.s, kind, datum, (NDDOSolutionUnrestricted) expSolution, true);
    }

    @Override
    protected void constructH() {
        h = new ParamHessianUnrestricted2((NDDOSolutionUnrestricted) opt.s, kind, datum, (NDDOSolutionUnrestricted) expSolution, true);
    }
}
