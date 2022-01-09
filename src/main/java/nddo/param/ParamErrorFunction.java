package nddo.param;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;

public abstract class ParamErrorFunction {
	protected final Solution soln;
	protected final NDDOAtom[] atoms;
	protected SimpleMatrix geomGradVector;
	protected Solution expSoln;
	protected NDDOAtom[] expAtoms;
	protected double HeatError, dipoleError, IEError, geomError, geomGradMag;

	public ParamErrorFunction(Solution soln, double refHeat) {
		this.atoms = soln.atoms;
		this.soln = soln;
		this.HeatError = (soln.hf - refHeat) * (soln.hf - refHeat);
		this.dipoleError = 0;
		this.IEError = 0;
	}

	public static ParamErrorFunction of(Solution s, double refHeat) {
		if (s instanceof SolutionR) return new ParamErrorFunctionR(s, refHeat);
		assert s instanceof SolutionU;
		return new ParamErrorFunctionU(s, refHeat);
	}

	public void createExpGeom(Solution expSoln) {
		this.expAtoms = expSoln.atoms;
		this.expSoln = expSoln;

		geomGradVector = new SimpleMatrix(expAtoms.length * 3, 1);
	}

	public void addDipoleError(double refDipole) {
		this.dipoleError = 400 * (soln.dipole - refDipole) * (soln.dipole - refDipole);
	}

	public void addIEError(double refIE) {
		this.IEError = 100 * (refIE + soln.homo) * (refIE + soln.homo);
	}

	// requires exp
	public void addGeomError() {
		double sum = 0;

		for (int i = 0, count = 0; i < expAtoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				double d = getGradient(i, j);
				geomGradVector.set(count, d);
				count++;
				sum += d * d;
			}
		}

		this.geomGradMag = Constants.KCAL * Math.sqrt(sum);
		this.geomError = 0.5 * Constants.KCAL * Constants.KCAL * sum;
	}

	public double getTotalError() {
		return HeatError + dipoleError + IEError + geomError;
	}

	public double getGeomGradMag() {
		return geomGradMag;
	}

	protected abstract double getGradient(int i, int j);
}
