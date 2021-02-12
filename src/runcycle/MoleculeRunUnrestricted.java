package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimizationUnrestricted;
import nddoparam.NDDOSolutionUnrestricted;
import nddoparam.param.NDDOParamGradientUnrestricted;
import nddoparam.param.NDDOParamHessianUnrestricted;

public class MoleculeRunUnrestricted extends MoleculeRun {
    public MoleculeRunUnrestricted(NDDOAtom[] atoms2, int charge, int mult, NDDOAtom[] expGeom, double[] datum, boolean runHessian, String trainingSet) {
        super(atoms2, charge, expGeom, datum, runHessian, trainingSet, mult);
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
        g = new NDDOParamGradientUnrestricted(atoms, charge, mult, Z, numit, (NDDOSolutionUnrestricted) opt.s);

    }

    protected void getH(int Z1, int param1, int Z2, int param2) {
        h = new NDDOParamHessianUnrestricted(atoms, charge, mult, Z1, param1, Z2, param2, (NDDOSolutionUnrestricted) opt.s);
    }

}
