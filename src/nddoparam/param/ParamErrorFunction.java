package nddoparam.param;

import nddoparam.NDDOAtom;
import nddoparam.Solution;
import nddoparam.SolutionR;
import nddoparam.SolutionU;
import scf.GTO;

import java.util.ArrayList;

public abstract class ParamErrorFunction {
	protected Solution soln, expSoln;
	protected double HeatError, dipoleError, IEError, geomError, geomGradient;
	protected NDDOAtom[] atoms, expAtoms;
	protected ArrayList<Double> bondErrors, angleErrors, bonds, angles,
			bondDerivatives,
			angleDerivatives;

	public ParamErrorFunction(Solution soln, double refHeat) {
		this.atoms = soln.atoms;
		this.soln = soln;
		this.HeatError = (soln.hf - refHeat) * (soln.hf - refHeat);
		this.dipoleError = 0;
		this.IEError = 0;

		this.angleErrors = new ArrayList<>();
		this.bondErrors = new ArrayList<>();
		this.bonds = new ArrayList<>();
		this.angles = new ArrayList<>();
		this.bondDerivatives = new ArrayList<>();
		this.angleDerivatives = new ArrayList<>();
	}

	public static ParamErrorFunction of(Solution s, double refHeat) {
		if (s instanceof SolutionR)
			return new ParamErrorFunctionR(s, refHeat);
		assert s instanceof SolutionU;
		return new ParamErrorFunctionU(s, refHeat);
	}

	private static double theta(double x1, double y1, double z1, double x2,
								double y2,
								double z2) {
		double sum = x1 * x2 + y1 * y2 + z1 * z2;
		double summ =
				mag(new double[]{x1, y1, z1}) * mag(new double[]{x2, y2, z2});
		return Math.acos(sum / summ) * 180 / Math.PI;
	}

	private static double[] cross(double[] a, double[] b) {
		return new double[]{a[1] * b[2] - a[2] * b[1],
				b[0] * a[2] - a[0] * b[2],
				a[0] * b[1] - b[0] * a[1]};
	}

	private static double[] normalizedVector(double[] v) {
		double val = mag(v);
		return new double[]{v[0] / val, v[1] / val, v[2] / val};
	}

	private static double mag(double[] vector) {
		double sum = 0;
		for (double v : vector) {
			sum += v * v;
		}
		return Math.sqrt(sum);
	}

	public void createExpGeom(Solution expSoln) {
		this.expAtoms = expSoln.atoms;
		this.expSoln = expSoln;
	}

	public void addDipoleError(double refDipole) {
		this.dipoleError =
				400 * (soln.dipole - refDipole) * (soln.dipole - refDipole);
	}

	public void addIEError(double refIE) {
		this.IEError = 100 * (refIE + soln.homo) * (refIE + soln.homo);
	}

	// requires exp
	public void addGeomError() {
		double sum = 0;
		for (int i = 0; i < expAtoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				double d = getGradient(i, j);
				sum += d * d;
			}
		}
		this.geomGradient = 627.5 * Math.sqrt(sum);
		this.geomError = 0.000049 * 627.5 * 627.5 * sum;
	}

	// requires exp
	@Deprecated
	public void addBondError(int atom1, int atom2, double ref) {
		double deriv = getDeriv(getR(atom1, atom2), atom2);
		this.bondDerivatives.add(deriv);
		this.bondErrors.add(0.49 * deriv * deriv);
	}

	// requires exp
	@Deprecated
	public void addAngleError(int atom1, int atom2, int atom3, double ref) {
		double[] coeff = getCoeff(atom1, atom2, atom3);
		double deriv = getDeriv(coeff, atom2);
		this.angleDerivatives.add(deriv);
		this.angleErrors.add(0.49 * deriv * deriv);
	}

	public double getTotalError() {
		double sum = HeatError + dipoleError + IEError + geomError;
		for (Double d : angleErrors) {
			sum += d;
		}
		for (Double d : bondErrors) {
			sum += d;
		}
		return sum;
	}

	public double getGeomGradient() {
		return geomGradient;
	}

	protected abstract double getDeriv(double[] coeff, int atom2);

	protected abstract double getGradient(int i, int j);

	protected double[] getR(int atom1, int atom2) {
		double length =
				GTO.R(atoms[atom1].getCoordinates(),
						atoms[atom2].getCoordinates()) /
						1.88973;
		this.bonds.add(length);
		double[] R = new double[]{
				expAtoms[atom2].getCoordinates()[0] -
						expAtoms[atom1].getCoordinates()[0],
				expAtoms[atom2].getCoordinates()[1] -
						expAtoms[atom1].getCoordinates()[1],
				expAtoms[atom2].getCoordinates()[2] -
						expAtoms[atom1].getCoordinates()[2]};
		double dist =
				GTO.R(expAtoms[atom2].getCoordinates(),
						expAtoms[atom1].getCoordinates());
		for (int x = 0; x < R.length; x++) {
			R[x] /= dist;
		}
		return R;
	}

	protected double[] getCoeff(int atom1, int atom2, int atom3) {
		double[] vector1 = new double[]{
				expAtoms[atom1].getCoordinates()[0] -
						expAtoms[atom2].getCoordinates()[0],
				expAtoms[atom1].getCoordinates()[1] -
						expAtoms[atom2].getCoordinates()[1],
				expAtoms[atom1].getCoordinates()[2] -
						expAtoms[atom2].getCoordinates()[2]};

		double[] vector2 = new double[]{
				expAtoms[atom3].getCoordinates()[0] -
						expAtoms[atom2].getCoordinates()[0],
				expAtoms[atom3].getCoordinates()[1] -
						expAtoms[atom2].getCoordinates()[1],
				expAtoms[atom3].getCoordinates()[2] -
						expAtoms[atom2].getCoordinates()[2]};

		double angle = theta(vector1[0], vector1[1], vector1[2], vector2[0],
				vector2[1],
				vector2[2]);
		this.angles.add(angle);
		double[] perpendicularVector =
				normalizedVector(cross(vector1, vector2));
		return normalizedVector(cross(vector2, perpendicularVector));
	}
}
