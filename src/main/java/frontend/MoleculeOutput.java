package frontend;

import runcycle.IMoleculeResult;
import runcycle.structs.RunnableMolecule;

public class MoleculeOutput { // MoleculeResult adapter
	public final RunnableMolecule updatedMolecule;
	public final long time;
	public final double hf, dipole, ie, geomGradient, totalError;
	public final ParamGradientOutput gradient;
	public final double[][] hessian;

	private MoleculeOutput(RunnableMolecule updatedMolecule, long time, double hf, double dipole, double ie,
						   double geomGradient, double totalError, ParamGradientOutput gradient,
						   double[][] hessian) {
		this.updatedMolecule = updatedMolecule;
		this.time = time;
		this.hf = hf;
		this.dipole = dipole;
		this.ie = ie;
		this.geomGradient = geomGradient;
		this.totalError = totalError;
		this.gradient = gradient;
		this.hessian = hessian;
	}

	public static MoleculeOutput from(IMoleculeResult result) {
		ParamGradientOutput pgo = new ParamGradientOutput(result.getHFDerivs(),
				result.getDipoleDerivs(), result.getIEDerivs(), result.getGeomDerivs(), result.getTotalGradients());

		return new MoleculeOutput(result.getUpdatedRm(), result.getTime(), result.getHF(),
				result.getDipole(), result.getIE(), result.getGeomGradient(), result.getTotalError(),
				pgo, result.getHessian());
	}

	public static class ParamGradientOutput {
		public final double[][] hf, dipole, ie, geom, total;

		public ParamGradientOutput(double[][] hf, double[][] dipole, double[][] ie, double[][] geom,
								   double[][] total) {
			this.hf = hf;
			this.dipole = dipole;
			this.ie = ie;
			this.geom = geom;
			this.total = total;
		}
	}
}
