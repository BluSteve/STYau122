package runcycle;

import nddoparam.*;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.*;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.input.RawMolecule;

import java.util.Arrays;
import java.util.concurrent.ExecutionException;

public class MoleculeRun {
	protected int[] atomTypes;
	protected double[] datum;
	protected NDDOAtom[] atoms, expGeom;
	protected Solution S, expS;
	protected GeometryOptimization opt;
	protected ParamGradient g;
	protected ParamHessian h;
	protected boolean isRunHessian, isExpAvail, restricted;
	protected int charge, mult;
	protected RawMolecule rawMolecule;
	protected long time;

	public MoleculeRun(RawMolecule rm, NDDOParams[] mp, int[] atomTypes,
					   boolean isRunHessian) {
		// todo change for am1
		atoms = RawMolecule.toMNDOAtoms(rm.atoms, (MNDOParams[]) mp);
		expGeom = rm.expGeom != null ? RawMolecule.toMNDOAtoms(rm.expGeom,
				(MNDOParams[]) mp) : null;
		charge = rm.charge;
		mult = rm.mult;
		datum = rm.datum.clone();
		rawMolecule = rm;
		restricted = rm.restricted;
		this.atomTypes = atomTypes;
		this.isRunHessian = isRunHessian;
		isExpAvail = expGeom != null;

		System.err.println(rawMolecule.index + " " + rawMolecule.name + " " +
				"started");
		try {
			run();
		} catch (Exception e) {
			System.err.println(
					"ERROR! " + Arrays.toString(e.getStackTrace()) + " " + rawMolecule.index + " " +
							rawMolecule.name);
//			try {
//				FileWriter fw = new FileWriter("errored.txt", true);
//				fw.write(e.getClass() + " " +
//						Arrays.toString(e.getStackTrace()) + " " +
//						rawMolecule.index + " " +
//						rawMolecule.name + "\n");
//				fw.close();
//			} catch (IOException ioException) {
//				ioException.printStackTrace();
//			}
		}
		System.err.println(rawMolecule.index + " " + rawMolecule.name +
				" finished in " + time);
	}

	public boolean isExpAvail() {
		return isExpAvail;
	}

	public RawMolecule getRawMolecule() {
		return rawMolecule;
	}

	public long getTime() {
		return time;
	}

	public double[] getDatum() {
		return datum;
	}

	public ParamGradient getG() {
		return g;
	}

	public ParamHessian getH() {
		return h;
	}

	private void run() throws ExecutionException, InterruptedException {
		StopWatch sw = new StopWatch();
		sw.start();

		opt = restricted ?
				new GeometryOptimizationR(atoms, charge) :
				new GeometryOptimizationU(atoms, charge, mult);
		S = getOpt().s; // NOT a clone

		// updates geom coords
		for (int i = 0; i < getAtoms().length; i++) {
			rawMolecule.atoms[i].coords = getAtoms()[i].getCoordinates();
		}

		if (this.getExpGeom() != null)
			expS = restricted ?
					new SolutionR(expGeom, charge) :
					new SolutionU(expGeom, charge, mult);

		g = restricted ?
				new ParamGradientR((SolutionR) getOpt().s, datum,
						(SolutionR) getExpS(), true, atomTypes) :
				new ParamGradientU((SolutionU) getOpt().s, datum,
						(SolutionU) getExpS(), false, atomTypes);
		h = restricted ?
				new ParamHessianR((ParamGradientR) g, true, atomTypes) :
				new ParamHessianU((ParamGradientU) g, false, atomTypes);


		g.computeGradients(); // time intensive step ~ 100 ms
		if (isRunHessian) {
			h.computeHessian(); // time intensive step ~ 700-800 ms
		}

		sw.stop();
		time = sw.getTime();
	}

	public Solution getS() {
		return S;
	}

	public int getCharge() {
		return charge;
	}

	public int getMult() {
		return mult;
	}

	public GeometryOptimization getOpt() {
		return opt;
	}

	public Solution getExpS() {
		return expS;
	}

	public NDDOAtom[] getAtoms() {
		return atoms;
	}

	public NDDOAtom[] getExpGeom() {
		return expGeom;
	}
}