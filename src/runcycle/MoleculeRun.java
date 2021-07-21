package runcycle;

import nddoparam.*;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamErrorFunction;
import nddoparam.param.ParamGradient;
import nddoparam.param.ParamHessian;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.input.RawMolecule;
import runcycle.output.OutputHandler;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

public class MoleculeRun implements MoleculeResult {
	protected double[] datum;
	protected NDDOAtom[] atoms, expGeom;
	protected Solution s, sExp;
	protected GeometryOptimization opt;
	protected ParamGradient g;
	protected ParamHessian h;
	protected boolean isRunHessian, isExpAvail, restricted;
	protected int charge, mult;
	protected RawMolecule rm;
	protected long time;

	public MoleculeRun(RawMolecule rm, NDDOParams[] mp, boolean isRunHessian) {
		// todo change for am1
		atoms = RawMolecule.toMNDOAtoms(rm.atoms, (MNDOParams[]) mp);
		expGeom = rm.expGeom != null ? RawMolecule.toMNDOAtoms(rm.expGeom,
				(MNDOParams[]) mp) : null;
		charge = rm.charge;
		mult = rm.mult;
		datum = rm.datum.clone();
		restricted = rm.restricted;
		this.rm = rm;
		this.isRunHessian = isRunHessian;
		isExpAvail = expGeom != null;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof MoleculeRun)) return false;
		else return ((MoleculeRun) obj).rm.index == this.rm.index;
	}

	public void run() {
		try {
			System.err.println(rm.index + " " + rm.name + " started");
			StopWatch sw = new StopWatch();
			sw.start();

			s = restricted ? new SolutionR(atoms, charge).setRm(rm) :
					new SolutionU(atoms, charge, mult).setRm(rm);

			opt = GeometryOptimization.of(s).compute();
			s = opt.getS();

			// updates geom coords
			for (int i = 0; i < atoms.length; i++) {
				rm.atoms[i].coords = atoms[i].getCoordinates();
			}

			if (expGeom != null) {
				sExp = restricted ?
						(new SolutionR(expGeom, charge)).setRm(rm) :
						(new SolutionU(expGeom, charge, mult)).setRm(rm);
			}

			g = ParamGradient.of(s, datum, sExp).compute();
			if (isRunHessian) h = ParamHessian.from(g).compute();

			sw.stop();
			time = sw.getTime();
			OutputHandler.outputOne(OutputHandler.toMoleculeOutput(this,
					isRunHessian),
					"dynamic-output");
			System.err.println(rm.index + " " + rm.name +
					" finished in " + time);
		} catch (Exception e) {
			e.printStackTrace();
			StringWriter errors = new StringWriter();
			e.printStackTrace(new PrintWriter(errors));

			String errorMessage = "ERROR! " + e.getClass() + " " +
					errors + " " + rm.index + " " + rm.name;
			logError(errorMessage);
		}
	}

	/**
	 * Logs the error and prevents this molecule from being run in the future
	 * by changing the isUsing parameter.
	 *
	 * @param errorMessage What to print to the console and log file.
	 */
	private void logError(String errorMessage) {
		System.err.println(errorMessage);
		rm.isUsing = false;
		try {
			FileWriter fw = new FileWriter("errored-molecules.txt", true);
			fw.write(errorMessage + "\n");
			fw.close();
		} catch (IOException ioException) {
			ioException.printStackTrace();
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
				"Hessian not found for molecule: " + rm.index + " " + rm.name);
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