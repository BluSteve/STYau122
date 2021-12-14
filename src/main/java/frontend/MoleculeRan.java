package frontend;

import runcycle.IMoleculeResult;
import runcycle.structs.RunnableMolecule;

public class MoleculeRan implements IMoleculeResult {
	MoleculeOutput mo;

	public MoleculeRan(MoleculeOutput mo) {
		this.mo = mo;
	}

	@Override
	public RunnableMolecule getUpdatedRm() {
		return mo.updatedMolecule;
	}

	@Override
	public long getTime() {
		return mo.time;
	}

	@Override
	public boolean isExpAvail() {
		return mo.updatedMolecule.expGeom != null;
	}

	@Override
	public double[][] getHessian() {
		if (mo.hessian != null) return mo.hessian;
		else throw new IllegalStateException(
				"Hessian not found for previously ran molecule: " + mo.updatedMolecule.debugName());
	}

	@Override
	public double getHF() {
		return mo.hf;
	}

	@Override
	public double getDipole() {
		return mo.dipole;
	}

	@Override
	public double getIE() {
		return mo.ie;
	}

	@Override
	public double getGeomGradient() {
		return mo.geomGradient;
	}

	@Override
	public double getTotalError() {
		return mo.totalError;
	}

	@Override
	public double[][] getHFDerivs() {
		return mo.gradient.hf;
	}

	@Override
	public double[][] getDipoleDerivs() {
		return mo.gradient.dipole;
	}

	@Override
	public double[][] getIEDerivs() {
		return mo.gradient.ie;
	}

	@Override
	public double[][] getGeomDerivs() {
		return mo.gradient.geom;
	}

	@Override
	public double[][] getTotalGradients() {
		return mo.gradient.total;
	}
}
