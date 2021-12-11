package runcycle;

import nddo.NDDOAtom;
import nddo.geometry.GeometryOptimization;
import nddo.solution.Solution;
import nddo.param.ParamErrorFunction;
import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.input.RawMolecule;
import runcycle.output.OutputHandler;

public class MoleculeRun implements MoleculeResult {
	protected double[] datum;
	protected NDDOAtom[] atoms, expGeom;
	protected Solution s, sExp;
	protected ParamGradient g;
	protected ParamHessian h;
	protected boolean isRunHessian, isExpAvail, restricted;
	protected int charge, mult;
	protected RawMolecule rm;
	protected long time;

	public MoleculeRun(RawMolecule rm, NDDOAtom[] atoms, NDDOAtom[] expGeom,
					   boolean isRunHessian) {
		// todo change for am1
		this.rm = rm;
		this.atoms = atoms;
		this.expGeom = expGeom;
		this.isRunHessian = isRunHessian;
		charge = rm.charge;
		mult = rm.mult;
		datum = rm.datum.clone();
		restricted = rm.restricted;
		isExpAvail = expGeom != null;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof MoleculeRun)) return false;
		else return ((MoleculeRun) obj).rm.index == this.rm.index;
	}

	public void run() {
		try {
			rm.getLogger().info("Started");
			StopWatch sw = new StopWatch();
			sw.start();


			s = GeometryOptimization.of(Solution.of(rm, atoms)).compute().getS();
			rm.getLogger().debug("Finished geometry optimization");

			if (isExpAvail) {
				sExp = Solution.of(rm, expGeom);
			}

			g = ParamGradient.of(s, datum, sExp).compute();
			rm.getLogger().debug("Finished param gradient");
			if (isRunHessian) h = ParamHessian.from(g).compute();
			rm.getLogger().debug("Finished param hessian");

			// updates geom coords, no IO involved
			for (int i = 0; i < s.atoms.length; i++) {
				rm.atoms[i].coords = s.atoms[i].getCoordinates();
			}


			sw.stop();
			time = sw.getTime();

			OutputHandler.outputOne(OutputHandler.toMoleculeOutput(this,
					isRunHessian), "dynamic-output");
			rm.getLogger().info("Finished in {}", time);
		} catch (Exception e) {
			rm.getLogger().error("", e);
		}
	}

	public boolean isExpAvail() {
		return isExpAvail;
	}

	public RawMolecule getRm() {
		return rm;
	}

	public long getTime() {
		return time;
	}

	public double[] getDatum() {
		return datum;
	}

	@Override
	public double[][] getHessian() {
		if (isRunHessian) return h.getHessian();
		else throw new IllegalStateException(
				"Hessian not found for molecule: " + rm.debugName());
	}

	@Override
	public double getHF() {
		return getS().hf;
	}

	@Override
	public double getDipole() {
		return getS().dipole;
	}

	@Override
	public double getIE() {
		return -getS().homo;
	}

	@Override
	public double getGeomGradient() {
		return getE().getGeomGradient();
	}

	@Override
	public double getTotalError() {
		return getE().getTotalError();
	}

	@Override
	public double[][] getHFDerivs() {
		return getG().getHFDerivs();
	}

	@Override
	public double[][] getDipoleDerivs() {
		return getG().getDipoleDerivs();
	}

	@Override
	public double[][] getIEDerivs() {
		return getG().getIEDerivs();
	}

	@Override
	public double[][] getGeomDerivs() {
		return getG().getGeomDerivs();
	}

	@Override
	public double[][] getTotalGradients() {
		return getG().getTotalGradients();
	}

	public ParamGradient getG() {
		return g;
	}

	public ParamHessian getH() {
		return h;
	}

	public Solution getS() {
		return s;
	}

	public ParamErrorFunction getE() {
		return g.getE();
	}
}