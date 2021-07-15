package runcycle;

import nddoparam.GeometryOptimization;
import nddoparam.NDDOAtom;
import nddoparam.Solution;
import nddoparam.param.ParamGradient;
import nddoparam.param.ParamHessian;
import runcycle.input.RawMolecule;

public abstract class MoleculeRun {
	protected int[] atomTypes;
	// TODO CHANGE THE ONES BELOW BACK TO PROTECTED
	protected NDDOAtom[] atoms;
	protected NDDOAtom[] expGeom;
	protected boolean isRunHessian;
	protected int charge;
	protected int size;
	protected int mult;
	protected GeometryOptimization opt;
	protected Solution expS;
	protected String kind;
	protected Solution S;
	protected RawMolecule rawMolecule;
	protected long time;
	protected double[] datum;
	protected ParamGradient g;
	protected ParamHessian h;
	protected boolean isExpAvail;

	public MoleculeRun(NDDOAtom[] atoms, int charge, NDDOAtom[] expGeom,
					   double[] datum,
					   boolean isRunHessian, String kind, int[] atomTypes,
					   int mult,
					   RawMolecule rawMolecule) {
		this.rawMolecule = rawMolecule;
		this.atomTypes = atomTypes;
		this.atoms = atoms;
		this.expGeom = expGeom;
		isExpAvail = expGeom != null;
		this.charge = charge;
		this.setRunHessian(isRunHessian);
		this.datum =
				datum; // reference heat + dipole + ionization. size = 1, 2
		// or 3
		this.size = atomTypes.length;
		this.mult = mult;
		this.setKind(kind);
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

	protected void routine() {
		runGradient(); // ~ 100 ms
		if (isRunHessian()) runHessian(); // ~ 700-800 ms
	}

	protected void generateGeomCoords() {
		for (int i = 0; i < getAtoms().length; i++) {
			rawMolecule.atoms[i].coords = getAtoms()[i].getCoordinates();
		}
	}

	protected void runHessian() {
		constructH();
		h.computeHessian(); // time intensive step
	}

	protected void runGradient() {
		constructG();
		g.computeGradients(); // time intensive step
	}

	protected abstract void constructG();

	public Solution getS() {
		return S;
	}

	protected abstract void constructH();

	public int getCharge() {
		return charge;
	}

	public int getSize() {
		return size;
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

	public boolean isRunHessian() {
		return isRunHessian;
	}

	public void setRunHessian(boolean runHessian) {
		isRunHessian = runHessian;
	}

	public String getKind() {
		return kind;
	}

	public void setKind(String kind) {
		this.kind = kind;
	}

	public NDDOAtom[] getAtoms() {
		return atoms;
	}

	public NDDOAtom[] getExpGeom() {
		return expGeom;
	}
}