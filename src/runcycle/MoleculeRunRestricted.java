package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimizationRestricted;
import nddoparam.NDDOParams;
import nddoparam.NDDOSolutionRestricted;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamGradientRestricted;
import nddoparam.param.ParamHessianRestricted;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.input.RawMolecule;

public class MoleculeRunRestricted extends MoleculeRun {
	public MoleculeRunRestricted(NDDOAtom[] atoms2, int charge, NDDOAtom[] expGeom,
								 double[] datum, boolean isRunHessian, String kind,
								 int[] atomTypes) {
		super(atoms2, charge, expGeom, datum, isRunHessian, kind, atomTypes, 1, null);
		metaRoutine();
	}

	public MoleculeRunRestricted(RawMolecule rm, NDDOParams[] mp, int[] atomTypes,
								 boolean isRunHessian) {
		super(RawMolecule.toMNDOAtoms(rm.atoms, (MNDOParams[]) mp), rm.charge,
				RawMolecule.toMNDOAtoms(rm.expGeom, (MNDOParams[]) mp),
				rm.datum, isRunHessian, rm.kind, atomTypes, 1, rm);
		metaRoutine();
	}

	private void metaRoutine() {
		StopWatch sw = new StopWatch();
		sw.start();

		opt = new NDDOGeometryOptimizationRestricted(atoms,
				charge); // ~ 74 ms for CH4 type d with expsoln

		generateGeomCoords();

		if (expGeom != null)
			expSolution = new NDDOSolutionRestricted(this.expGeom, charge); // ~ 60 ms

		routine(); // ~700-800 ms

		sw.stop();
		time = sw.getTime();
	}

	@Override
	protected void constructG() {
		g = new ParamGradientRestricted((NDDOSolutionRestricted) opt.s, kind, datum,
				(NDDOSolutionRestricted) expSolution, true);
	}

	@Override
	protected void constructH() {
		h = new ParamHessianRestricted((ParamGradientRestricted) g, true);
	}
}
