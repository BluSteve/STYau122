package runcycle;

import nddoparam.*;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.*;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.input.RawMolecule;
import runcycle.output.OutputHandler;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.concurrent.*;

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

	public void run() {
		try {
			System.err
					.println(rawMolecule.index + " " + rawMolecule.name + " " +
							"started");
			ExecutorService executorService =
					Executors.newSingleThreadExecutor();
			Future future = executorService.submit(() -> {
				StopWatch sw = new StopWatch();
				sw.start();

				opt = restricted ?
						new GeometryOptimizationR(atoms, charge) :
						new GeometryOptimizationU(atoms, charge, mult);
				S = getOpt().s; // NOT a clone

				// updates geom coords
				for (int i = 0; i < getAtoms().length; i++) {
					rawMolecule.atoms[i].coords =
							getAtoms()[i].getCoordinates();
				}

				if (this.getExpGeom() != null)
					expS = restricted ?
							new SolutionR(expGeom, charge) :
							new SolutionU(expGeom, charge, mult);

				g = restricted ?
						new ParamGradientR((SolutionR) S, datum,
								(SolutionR) expS, true) :
						new ParamGradientU((SolutionU) S, datum,
								(SolutionU) expS, false);
				h = restricted ?
						new ParamHessianR((ParamGradientR) g, true,
								atomTypes) :
						new ParamHessianU((ParamGradientU) g, false,
								atomTypes);


				g.compute(); // time intensive step ~ 100 ms
				if (isRunHessian) {
					h.compute(); // time intensive step ~ 700-800 ms
				}

				sw.stop();
				time = sw.getTime();
				OutputHandler.outputOne(OutputHandler.toMoleculeOutput(this),
						"dynamic-output");
				System.err.println(rawMolecule.index + " " + rawMolecule.name +
						" finished in " + time);
			});

			try {
				future.get(600, TimeUnit.SECONDS);
			} catch (TimeoutException e) {
				future.cancel(true);

				System.err.println("TIMEOUT! " + rawMolecule.index + " " +
						rawMolecule.name);

				rawMolecule.isUsing = false;

				try {
					FileWriter fw = new FileWriter("errored-molecules.txt",
							true);
					fw.write("TIMEOUT! " + rawMolecule.index + " " +
							rawMolecule.name + "\n");
					fw.close();
				} catch (IOException ioException) {
					ioException.printStackTrace();
				}
			} finally {
				executorService.shutdown();
			}
		} catch (Exception e) {
			System.err.println(
					"ERROR! " + e.getClass() + " " +
							Arrays.toString(e.getStackTrace()) + " " +
							rawMolecule.index + " " +
							rawMolecule.name);
			rawMolecule.isUsing = false;
			try {
				FileWriter fw = new FileWriter("errored-molecules.txt", true);
				fw.write("ERROR! " + e.getClass() + " " +
						Arrays.toString(e.getStackTrace()) + " " +
						rawMolecule.index + " " +
						rawMolecule.name + "\n");
				fw.close();
			} catch (IOException ioException) {
				ioException.printStackTrace();
			}
		}
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