package frontend.json;

import runcycle.IMoleculeResult;
import runcycle.structs.RunnableMolecule;

public class MoleculeSerial implements IMoleculeResult {
	public RunnableMolecule updatedRm;
	public long time;
	public double hf, dipole, ie, geomGradient, totalError;
	public double[][] hfpg, dipolepg, iepg, geompg, totalpg, hessian;

	public MoleculeSerial() {
	}

	private MoleculeSerial(RunnableMolecule updatedRm, long time, double hf, double dipole, double ie,
						   double geomGradient, double totalError, double[][] hfpg, double[][] dipolepg,
						   double[][] iepg, double[][] geompg, double[][] totalpg, double[][] hessian) {
		this.updatedRm = updatedRm;
		this.time = time;
		this.hf = hf;
		this.dipole = dipole;
		this.ie = ie;
		this.geomGradient = geomGradient;
		this.totalError = totalError;
		this.hfpg = hfpg;
		this.dipolepg = dipolepg;
		this.iepg = iepg;
		this.geompg = geompg;
		this.totalpg = totalpg;
		this.hessian = hessian;
	}

	public static MoleculeSerial from(IMoleculeResult result) {
		return new MoleculeSerial(result.getUpdatedRm(), result.getTime(), result.getHF(),
				result.getDipole(), result.getIE(), result.getGeomGradient(), result.getTotalError(),
				result.getHFDerivs(), result.getDipoleDerivs(), result.getIEDerivs(), result.getGeomDerivs(),
				result.getTotalGradients(), result.getHessian());
	}

	@Override
	public RunnableMolecule getUpdatedRm() {
		return updatedRm;
	}

	@Override
	public long getTime() {
		return time;
	}

	@Override
	public boolean isExpAvail() {
		return updatedRm.expGeom != null;
	}

	@Override
	public double[][] getHessian() {
		if (hessian != null) return hessian;
		else throw new IllegalStateException(
				"Hessian not found for previously ran molecule: " + updatedRm.debugName());
	}

	@Override
	public double getHF() {
		return hf;
	}

	@Override
	public double getDipole() {
		return dipole;
	}

	@Override
	public double getIE() {
		return ie;
	}

	@Override
	public double getGeomGradient() {
		return geomGradient;
	}

	@Override
	public double getTotalError() {
		return totalError;
	}

	@Override
	public double[][] getHFDerivs() {
		return hfpg;
	}

	@Override
	public double[][] getDipoleDerivs() {
		return dipolepg;
	}

	@Override
	public double[][] getIEDerivs() {
		return iepg;
	}

	@Override
	public double[][] getGeomDerivs() {
		return geompg;
	}

	@Override
	public double[][] getTotalGradients() {
		return totalpg;
	}
}
