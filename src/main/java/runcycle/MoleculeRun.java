package runcycle;

import nddo.NDDOAtom;
import nddo.geometry.GeometryOptimization;
import nddo.param.ParamErrorFunction;
import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import nddo.solution.Solution;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.structs.Atom;
import runcycle.structs.RunnableMolecule;

public final class MoleculeRun implements IMoleculeResult {
	private final NDDOAtom[] nddoAtoms, expGeom;
	private final boolean withHessian, isExpAvail;
	private final RunnableMolecule rm;
	private final double[] datum;
	private Solution s, sExp;
	private ParamGradient g;
	private ParamHessian h;
	private Atom[] newAtoms;
	private long time;

	public MoleculeRun(RunnableMolecule rm, NDDOAtom[] nddoAtoms, NDDOAtom[] expGeom, double[] datum,
					   boolean withHessian) {
		this.rm = rm;
		this.nddoAtoms = nddoAtoms;
		this.expGeom = expGeom;
		this.datum = datum;
		this.withHessian = withHessian;

		isExpAvail = expGeom != null;
	}

	public void run() {
		try {
			rm.getLogger().info("Started");
			StopWatch sw = new StopWatch();
			sw.start();


			s = GeometryOptimization.of(Solution.of(rm, nddoAtoms)).compute().getS();
			rm.getLogger().debug("Finished geometry optimization");

			if (isExpAvail) {
				sExp = Solution.of(rm, expGeom);
			}

			g = ParamGradient.of(s, datum, sExp).compute();
			rm.getLogger().debug("Finished param gradient");
			if (withHessian) h = ParamHessian.from(g).compute();
			rm.getLogger().debug("Finished param hessian");

			// stores new optimized geometry
			newAtoms = new Atom[s.atoms.length];
			for (int i = 0; i < newAtoms.length; i++) {
				newAtoms[i] = new Atom(s.atoms[i].getAtomProperties().getZ(), s.atoms[i].getCoordinates());
			}


			sw.stop();
			time = sw.getTime();

			rm.getLogger().info("Finished in {}", time);
		} catch (Exception e) {
			rm.getLogger().error("", e);
		}
	}

	public boolean isExpAvail() {
		return isExpAvail;
	}

	public RunnableMolecule getUpdatedRm() {
		return new RunnableMolecule(rm, newAtoms, rm.expGeom, rm.datum);
	}

	public long getTime() {
		return time;
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

	@Override
	public double[][] getHessian() {
		if (withHessian) return h.getHessian();
		else throw new IllegalStateException("Hessian not found for molecule: " + rm.debugName());
	}

	public ParamGradient getG() {
		return g;
	}

	public ParamHessian getH() {
		return h;
	}

	/**
	 * @return Original, unoptimized Solution object.
	 */
	public Solution getS() {
		return s;
	}

	public ParamErrorFunction getE() {
		return g.getE();
	}
}