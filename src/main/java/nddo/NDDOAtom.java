package nddo;

import nddo.geometry.GeometryDerivative;
import nddo.geometry.GeometrySecondDerivative;
import nddo.param.ParamDerivative;
import scf.Atom;
import scf.AtomProperties;
import scf.OrbitalProperties;

public abstract class NDDOAtom extends Atom {
	protected static final double bohr = 1.88973;
	protected final NDDOParams np;
	public double p0, p1, p2, D1, D2;
	protected NDDO6G[] orbitals;

	public NDDOAtom(AtomProperties atomProperties, double[] coordinates, NDDOParams np) {
		super(atomProperties, coordinates);
		this.np = np.clone();
		this.p0 = p0();
		this.D1 = D1();
		this.D2 = D2();
		this.p1 = p1();
		this.p2 = p2();

		if (this.atomProperties.getZ() == 1) {
			this.p1 = 0;
			this.p2 = 0;
			this.D1 = 0;
			this.D2 = 0;
		}

		OrbitalProperties[] orbitalProperties = this.atomProperties.getOrbitals();
		orbitals = new NDDO6G[orbitalProperties.length];
		for (int x = 0; x < orbitals.length; x++) {
			switch (orbitalProperties[x].getType()) {
				case "s":
					orbitals[x] = new NDDO6G(this, orbitalProperties[x], np.getZetas(), np.getBetas(), np.getUss());
					break;
				case "p":
					orbitals[x] = new NDDO6G(this, orbitalProperties[x], np.getZetap(), np.getBetap(), np.getUpp());
					break;
			}
		}
	}

	/**
	 * Returns a brand new NDDOAtom object. Everything is pass-by-value.
	 * @param np New NDDOParams.
	 * @return New NDDOAtom.
	 */
	public abstract NDDOAtom withNewParams(NDDOParams np);

	/**
	 * Returns a brand new NDDOAtom object. Everything is pass-by-value.
	 * @param coordinates New coordinates.
	 * @return New NDDOAtom.
	 */
	public abstract NDDOAtom withNewCoords(double[] coordinates);

	public NDDO6G[] getOrbitals() {
		return this.orbitals;
	}

	public double getMass() {
		return atomProperties.getMass();
	}

	public double getHeat() {
		return atomProperties.getHeat();
	}

	/**
	 * @return Cloned NDDOParams.
	 */
	public NDDOParams getParams() {
		return np.clone();
	}

	public NDDO6G s() {
		return this.orbitals[0]; // todo what is this
	}

	public double V(NDDO6G a, NDDO6G b) {
		return -this.atomProperties.getQ() * NDDO6G.getG(a, b, this.s(), this.s());
	}

	public double Vderiv(NDDO6G a, NDDO6G b, int tau) {
		return -this.atomProperties.getQ() * GeometryDerivative.getGderiv(a, b, this.s(), this.s(), tau);
	}

	public double Vderiv2(NDDO6G a, NDDO6G b, int tau1, int tau2) {
		return -this.atomProperties.getQ() * GeometrySecondDerivative.getGderiv2(a, b, this.s(), this.s(), tau1, tau2);
	}

	public double VParamDeriv(NDDO6G a, NDDO6G b, int num, int type) {
		return -this.atomProperties.getQ() * ParamDerivative.getGderiv(this.s(), this.s(), a, b, num, type);
	}

	protected double p0() {
		return 27.2114 / (2 * np.getGss());
	}

	protected double D1() {
		return (2 * atomProperties.getPeriod() + 1) / Math.sqrt(3) *
				Math.pow(4 * np.getZetas() * np.getZetap(), atomProperties.getPeriod() + 0.5) /
				Math.pow(np.getZetas() + np.getZetap(), 2 * atomProperties.getPeriod() + 2);
	}

	protected double D2() {
		return 1 / np.getZetap() *
				Math.sqrt((2 * atomProperties.getPeriod() + 1) * (2 * atomProperties.getPeriod() + 2) / 20.0);
	}

	protected double p1() {
		double guess = 0;
		double newguess = 0.5 * Math.pow(D1 * D1 * 27.2114 / (np.getHsp()), 1.0 / 3);
		while (Math.abs(guess - newguess) > 1E-12) {
			guess = newguess;
			double f = 1 / guess - 1 / Math.sqrt(guess * guess + D1 * D1) - 4 * np.getHsp() / 27.2114;
			double fprime = -1 / (guess * guess) + guess / Math.pow(guess * guess + D1 * D1, 1.5);
			newguess = guess - f / fprime;
		}
		return newguess;
	}

	protected double p2() {
		double guess = 0;
		double newguess = 0.5;
		while (Math.abs(guess - newguess) > 1E-12) {
			guess = newguess;
			double f = 1 / guess + 1 / Math.sqrt(guess * guess + 2 * D2 * D2) -
							2 / Math.sqrt(guess * guess + D2 * D2) - 4 * (np.getGpp() - np.getGp2()) / 27.2114;
			double fprime = -1 / (guess * guess) - guess / Math.pow(guess * guess + 2 * D2 * D2, 1.5) +
					2 * guess / Math.pow(guess * guess + D2 * D2, 1.5);
			newguess = guess - f / fprime;
		}
		return newguess;
	}

	public abstract double crf(NDDOAtom b);

	public abstract double crfDeriv(NDDOAtom b, int tau);

	public abstract double crfDeriv2(NDDOAtom b, int tau1, int tau2);

	public abstract double crfAlphaDeriv(NDDOAtom b, int num);

	@Override
	public abstract NDDOAtom clone();
}
