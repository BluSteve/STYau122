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
import java.util.Arrays;
import java.util.concurrent.*;

public class MoleculeRun {
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

	public void run() {
		try {
			System.err.println(rm.index + " " + rm.name + " started");
			ExecutorService executorService =
					Executors.newSingleThreadExecutor();
			Future<?> future = executorService.submit(() -> {
				StopWatch sw = new StopWatch();
				sw.start();

				opt = restricted ?
						new GeometryOptimizationR(atoms, charge) :
						new GeometryOptimizationU(atoms, charge, mult);
				s = opt.s.setRm(rm); // NOT a clone

				// updates geom coords
				for (int i = 0; i < atoms.length; i++) {
					rm.atoms[i].coords = atoms[i].getCoordinates();
				}

				if (expGeom != null)
					sExp = restricted ?
							(new SolutionR(expGeom, charge)).setRm(rm) :
							(new SolutionU(expGeom, charge, mult)).setRm(rm);

				g = ParamGradient.of(s, datum, sExp).compute();
				if (isRunHessian) h = ParamHessian.from(g).compute();

				sw.stop();
				time = sw.getTime();
				OutputHandler.outputOne(OutputHandler.toMoleculeOutput(this),
						"dynamic-output");
				System.err.println(rm.index + " " + rm.name +
						" finished in " + time);
			});

			try {
				future.get(600, TimeUnit.SECONDS);
			} catch (TimeoutException e) {
				e.printStackTrace();
				future.cancel(true);

				String timeoutMessage = "TIMEOUT! " + rm.index + " " + rm.name;

				logError(timeoutMessage);
			} finally {
				executorService.shutdown();
			}
		} catch (Exception e) {
			e.printStackTrace();
			String errorMessage = "ERROR! " + e.getClass() + " " +
					Arrays.toString(e.getStackTrace()) + " " +
					rm.index + " " + rm.name;

			logError(errorMessage);
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

	public ParamGradient getG() {
		return g;
	}

	public ParamHessian getH() {
		return h;
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

	public Solution getS() {
		return s;
	}

	public ParamErrorFunction getE() {
		return g.getE();
	}
}