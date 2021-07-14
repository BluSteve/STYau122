package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimizationUnrestricted;
import nddoparam.NDDOParams;
import nddoparam.NDDOSolutionUnrestricted;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamGradientUnrestricted2;
import nddoparam.param.ParamHessianUnrestricted2;
import runcycle.input.RawMolecule;

public class MoleculeRunUnrestricted extends MoleculeRun {
    public MoleculeRunUnrestricted(NDDOAtom[] atoms2, int charge, int mult, NDDOAtom[] expGeom, double[] datum, boolean isRunHessian, String kind, int[] atomTypes) {
        super(atoms2, charge, expGeom, datum, isRunHessian, kind, atomTypes, mult);

        opt = new NDDOGeometryOptimizationUnrestricted(atoms, charge, mult);

        newGeomCoords = "UHF\n" + "CHARGE=" + charge + "\nMULT=" + mult + "\n";
        generateGeomCoords();
        if (this.expGeom != null) expSolution = new NDDOSolutionUnrestricted(this.expGeom, charge, mult);

        routine();
    }

    public MoleculeRunUnrestricted(RawMolecule rm, NDDOParams[] mp, int[] atomTypes, boolean isRunHessian) {
        this(RawMolecule.toMNDOAtoms(rm.atoms, (MNDOParams[]) mp), rm.charge, rm.mult, RawMolecule.toMNDOAtoms(rm.expGeom, (MNDOParams[]) mp),
                rm.datum, isRunHessian, rm.kind, atomTypes);
    }

    @Override
    protected void constructG() {
        // TODO Change this line's hardcoding once analytical has been implemented for unrestricted.
        g = new ParamGradientUnrestricted2((NDDOSolutionUnrestricted) opt.s, kind, datum, (NDDOSolutionUnrestricted) expSolution, false);
    }

    @Override
    // TODO Change this line's hardcoding once analytical has been implemented for unrestricted.
    protected void constructH() {
        h = new ParamHessianUnrestricted2((ParamGradientUnrestricted2) g, false);
    }
}
