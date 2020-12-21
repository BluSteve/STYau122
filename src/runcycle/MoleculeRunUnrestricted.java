package runcycle;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDOGeometryOptimizationUnrestricted;
import mndoparam.mndo.MNDOSolutionUnrestricted;
import mndoparam.param.MNDOParamGradientUnrestricted;
import mndoparam.param.MNDOParamHessianUnrestricted;

public class MoleculeRunUnrestricted extends MoleculeRun {
    public MoleculeRunUnrestricted(MNDOAtom[] atoms2, int charge, int mult, MNDOAtom[] expGeom, double[] datum, boolean runHessian, String trainingSet) {
        super(atoms2, charge, expGeom, datum, runHessian, trainingSet, mult);
        opt = new MNDOGeometryOptimizationUnrestricted(atoms, charge, mult);

        newGeomCoords = "UHF\n" + "CHARGE=" + charge + "\nMULT=" + mult + "\n";
        generateGeomCoords();
        if (this.expGeom != null) expSolution = new MNDOSolutionUnrestricted(this.expGeom, charge, mult);

        if (runHessian) {
            runHessian();
        } else {
            hessianStr = "";
        }
        runGradient();

        outputErrorFunction();
    }

    protected void getG(int Z, int numit, int mult) {
        g = new MNDOParamGradientUnrestricted(atoms, charge, mult, Z, numit, (MNDOSolutionUnrestricted) opt.s);

    }

    protected void getH(int Z1, int param1, int Z2, int param2) {
        h = new MNDOParamHessianUnrestricted(atoms, charge, mult, Z1, param1, Z2, param2, (MNDOSolutionUnrestricted) opt.s);
    }

}
