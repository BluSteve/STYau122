package runcycle;

import nddoparam.GeometryOptimization;
import nddoparam.NDDOAtom;
import nddoparam.Solution;
import nddoparam.param.ParamGradient;
import nddoparam.param.ParamHessian;
import runcycle.input.RawMolecule;

public abstract class MoleculeRun {
	public String[] output;
	public String hessianStr, newGeomCoords;
	public String totalHeatDeriv = "";
	public String totalIonizationDeriv = "";
	public String totalDipoleDeriv = "";
	public String totalGeomDeriv = "";
	public String totalExcelStr = "";
	public String excelStr = "";
	public String excelStr2 = "";
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

	public boolean isExpAvail() {
		return isExpAvail;
	}

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

	private static void appendToSB(double[] array, StringBuilder sb) {
		if (sb.length() > 0) sb.append(',');
		for (int i = 0; i < array.length - 1; i++)
			sb.append(array[i]).append(',');
		sb.append(array[array.length - 1]);
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
		else hessianStr = "";

		outputErrorFunction();
	}

	protected void generateGeomCoords() {
		for (int i = 0; i < getAtoms().length; i++) {
			rawMolecule.atoms[i].coords = getAtoms()[i].getCoordinates();
		}
	}

	protected void runHessian() {
		constructH();
		h.computeHessian(); // time intensive step
		StringBuilder hessianSB = new StringBuilder();
		appendToSB(h.getHessianUT(), hessianSB);
		hessianStr = hessianSB.toString();
	}

	protected void runGradient() {
		constructG();
		g.computeGradients(); // time intensive step
		StringBuilder HFDerivsSB =
				new StringBuilder(datum[0] + "," + g.getS().hf);
		StringBuilder dipoleDerivsSB =
				new StringBuilder(datum[1] + "," + g.getS().dipole);
		StringBuilder IEDerivsSB =
				new StringBuilder(datum[2] + "," + -g.getS().homo);
		StringBuilder geomDerivsSB =
				new StringBuilder("0," + g.getE().geomGradient);
		StringBuilder mainDataSB = new StringBuilder();
		mainDataSB.append(g.getE().getTotalError());

		for (int i = 0; i < g.getS().getUniqueZs().length; i++) {
			switch (getKind()) {
				case "a":
					appendToSB(g.depad(g.getHFDerivs())[i], HFDerivsSB);
					break;
				case "b":
					appendToSB(g.depad(g.getHFDerivs())[i], HFDerivsSB);
					appendToSB(g.depad(g.getDipoleDerivs())[i],
							dipoleDerivsSB);
					break;
				case "c":
					appendToSB(g.depad(g.getHFDerivs())[i], HFDerivsSB);
					appendToSB(g.depad(g.getIEDerivs())[i], IEDerivsSB);
					break;
				case "d":
					appendToSB(g.depad(g.getHFDerivs())[i], HFDerivsSB);
					appendToSB(g.depad(g.getDipoleDerivs())[i],
							dipoleDerivsSB);
					appendToSB(g.depad(g.getIEDerivs())[i], IEDerivsSB);
					break;
			}
			if (getExpS() != null)
				appendToSB(g.depad(g.getGeomDerivs())[i], geomDerivsSB);
			appendToSB(g.depad(g.getTotalGradients())[i], mainDataSB);
		}

		appendToSB(
				new double[]{datum[0], g.getS().hf, datum[1], g.getS().dipole,
						datum[2], g.getS().homo, 0, g.getE().geomGradient},
				mainDataSB);

		totalHeatDeriv = HFDerivsSB.toString();
		totalDipoleDeriv = dipoleDerivsSB.toString();
		totalIonizationDeriv = IEDerivsSB.toString();
		totalGeomDeriv = geomDerivsSB.toString();
		totalExcelStr = mainDataSB.toString();
	}

	protected void outputErrorFunction() {
		System.out.println("Error function: " + g.getE().getTotalError());
		output = new String[]{totalExcelStr, totalHeatDeriv, totalDipoleDeriv,
				totalIonizationDeriv, totalGeomDeriv};
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