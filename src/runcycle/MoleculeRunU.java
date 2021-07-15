package runcycle;

import nddoparam.GeometryOptimizationU;
import nddoparam.NDDOAtom;
import nddoparam.NDDOParams;
import nddoparam.SolutionU;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamGradientU;
import nddoparam.param.ParamHessianU;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.input.RawMolecule;

public class MoleculeRunU extends MoleculeRun {
	public MoleculeRunU(NDDOAtom[] atoms2, int charge, int mult,
						NDDOAtom[] expGeom,
						double[] datum, boolean isRunHessian, String kind,
						int[] atomTypes) {
		super(atoms2, charge, expGeom, datum, isRunHessian, kind, atomTypes,
				mult, null);
		metaRoutine();
	}

	public MoleculeRunU(RawMolecule rm, NDDOParams[] mp, int[] atomTypes,
						boolean isRunHessian) {
		super(RawMolecule.toMNDOAtoms(rm.atoms, (MNDOParams[]) mp), rm.charge,
				RawMolecule.toMNDOAtoms(rm.expGeom, (MNDOParams[]) mp),
				rm.datum, isRunHessian, rm.kind, atomTypes, rm.mult, rm);
		metaRoutine();
	}

	private void metaRoutine() {
		StopWatch sw = new StopWatch();
		sw.start();

		opt = new GeometryOptimizationU(getAtoms(), getCharge(), getMult());

		generateGeomCoords();
		if (this.getExpGeom() != null)
			expS = new SolutionU(this.getExpGeom(), getCharge(), getMult());

		routine();

		sw.stop();
		time = sw.getTime();
	}

	@Override
	protected void constructG() {
		// TODO Change this line's hardcoding once analytical has been
		//  implemented for
		//  unrestricted.
		g = new ParamGradientU((SolutionU) getOpt().s, getKind(), datum,
				(SolutionU) getExpS(), false);
	}

	@Override
	// TODO Change this line's hardcoding once analytical has been implemented
	//  for
	//  unrestricted.
	protected void constructH() {
		h = new ParamHessianU((ParamGradientU) g, false);
	}
}
