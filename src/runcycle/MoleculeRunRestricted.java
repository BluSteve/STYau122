package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimizationRestricted;
import nddoparam.NDDOSolutionRestricted;
import nddoparam.param.NDDOParamGradientRestricted;
import nddoparam.param.NDDOParamHessianRestricted;

public class MoleculeRunRestricted extends MoleculeRun {
    public MoleculeRunRestricted(NDDOAtom[] atoms2, int charge, NDDOAtom[] expGeom, double[] datum, boolean runHessian, String trainingSet) {
        super(atoms2, charge, expGeom, datum, runHessian, trainingSet, 1);
        opt = new NDDOGeometryOptimizationRestricted(atoms, charge);

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

    protected void getG(int Z, int numit, int mult) { // TODO make everything use NDDOSolution instead
        g = new NDDOParamGradientRestricted(atoms, charge, Z, numit, (NDDOSolutionRestricted) opt.s);

    }

    protected void getH(int Z1, int param1, int Z2, int param2) {
        h = new NDDOParamHessianRestricted(atoms, charge, Z1, param1, Z2, param2, (NDDOSolutionRestricted) opt.s);
    }

}
