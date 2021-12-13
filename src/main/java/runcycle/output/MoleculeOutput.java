package runcycle.output;

import runcycle.MoleculeResult;
import runcycle.structs.RunnableMolecule;

public class MoleculeOutput {
	public final RunnableMolecule runnableMolecule;
	public final long time;
	public final double hf, dipole, ie, geomGradient, totalError;
	public final ParamGradientOutput gradient;
	public final double[][] hessian;

	private MoleculeOutput(RunnableMolecule runnableMolecule, long time, double hf, double dipole, double ie,
						   double geomGradient, double totalError, ParamGradientOutput gradient,
						   double[][] hessian) {
		this.runnableMolecule = runnableMolecule;
		this.time = time;
		this.hf = hf;
		this.dipole = dipole;
		this.ie = ie;
		this.geomGradient = geomGradient;
		this.totalError = totalError;
		this.gradient = gradient;
		this.hessian = hessian;
	}

	public static MoleculeOutput from(MoleculeResult result) {
		ParamGradientOutput pgo = new ParamGradientOutput(result.getHFDerivs(),
				result.getDipoleDerivs(), result.getIEDerivs(), result.getGeomDerivs(), result.getTotalGradients());

		return new MoleculeOutput(result.getRm(), result.getTime(), result.getHF(),
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
