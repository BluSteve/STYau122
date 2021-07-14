package runcycle;

import nddoparam.GeometryOptimizationR;
import nddoparam.NDDOAtom;
import nddoparam.NDDOParams;
import nddoparam.SolutionR;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamGradientR;
import nddoparam.param.ParamHessianR;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.input.RawMolecule;

public class MoleculeRunR extends MoleculeRun {
	public MoleculeRunR(NDDOAtom[] atoms2, int charge, NDDOAtom[] expGeom,
						double[] datum, boolean isRunHessian, String kind,
						int[] atomTypes) {
		super(atoms2, charge, expGeom, datum, isRunHessian, kind, atomTypes, 1, null);
		metaRoutine();
	}

	public MoleculeRunR(RawMolecule rm, NDDOParams[] mp, int[] atomTypes,
						boolean isRunHessian) {
		super(RawMolecule.toMNDOAtoms(rm.atoms, (MNDOParams[]) mp), rm.charge,
				RawMolecule.toMNDOAtoms(rm.expGeom, (MNDOParams[]) mp),
				rm.datum, isRunHessian, rm.kind, atomTypes, 1, rm);
		metaRoutine();
	}

	private void metaRoutine() {
		StopWatch sw = new StopWatch();
		sw.start();

		opt = new GeometryOptimizationR(atoms,
				charge); // ~ 74 ms for CH4 type d with expsoln

		generateGeomCoords();

		if (expGeom != null)
			expSolution = new SolutionR(this.expGeom, charge); // ~ 60 ms

		routine(); // ~700-800 ms

		sw.stop();
		time = sw.getTime();
	}

	@Override
	protected void constructG() {
		g = new ParamGradientR((SolutionR) opt.s, kind, datum,
				(SolutionR) expSolution, true);
	}

	@Override
	protected void constructH() {
		h = new ParamHessianR((ParamGradientR) g, true);
	}
}
