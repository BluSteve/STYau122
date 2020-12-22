package runcycle;

import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOGeometryOptimizationRestricted;
import nddoparam.mndo.MNDOSolutionRestricted;
import nddoparam.param.MNDOParamGradientRestricted;
import nddoparam.param.MNDOParamHessianRestricted;

public class MoleculeRunRestricted extends MoleculeRun {
    public MoleculeRunRestricted(MNDOAtom[] atoms2, int charge, MNDOAtom[] expGeom, double[] datum, boolean runHessian, String trainingSet) {
        super(atoms2, charge, expGeom, datum, runHessian, trainingSet, 1);
        opt = new MNDOGeometryOptimizationRestricted(atoms, charge);

        super.newGeomCoords = "RHF\n" + "CHARGE=" + charge + "\nMULT=1\n";
        super.generateGeomCoords();
        if (expGeom != null) expSolution = new MNDOSolutionRestricted(this.expGeom, charge);

        if (runHessian) {
            runHessian();
        } else {
            hessianStr = "";
        }
        runGradient();

        outputErrorFunction();
    }

    protected void getG(int Z, int numit, int mult) { // TODO make everything use MNDOSolution instead
        g = new MNDOParamGradientRestricted(atoms, charge, Z, numit, (MNDOSolutionRestricted) opt.s);

    }

    protected void getH(int Z1, int param1, int Z2, int param2) {
        h = new MNDOParamHessianRestricted(atoms, charge, Z1, param1, Z2, param2, (MNDOSolutionRestricted) opt.s);
    }

}
