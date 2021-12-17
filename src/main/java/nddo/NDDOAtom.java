package nddo;

import nddo.structs.AtomProperties;

public abstract class NDDOAtom {
	protected final double[] coordinates;
	protected final AtomProperties atomProperties;
	protected final NDDOParams np;
	public double p0, p1, p2, D1, D2;
	protected NDDOOrbital[] orbitals;

	/**
	 * The standard representation of an atom used all over nddo.
	 *
	 * @param atomProperties Fixed atom properties.
	 * @param coordinates    Coordinates will be cloned/ are passed-by-value.
	 * @param np             NDDOParams will be cloned/ are pass-by-value.
	 */
	public NDDOAtom(AtomProperties atomProperties, double[] coordinates, NDDOParams np) {
		this.atomProperties = atomProperties;
		this.coordinates = coordinates.clone();
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
	}

	/**
	 * Get coordinates object with no cloning. Can be modified in place!
	 *
	 * @return Original coordinates array object.
	 */
	public double[] getCoordinates() {
		return this.coordinates;
	}

	public AtomProperties getAtomProperties() {
		return this.atomProperties;
	}

	/**
	 * Returns a brand new NDDOAtom object. Everything is pass-by-value.
	 *
	 * @param np New NDDOParams.
	 * @return New NDDOAtom.
	 */
	public abstract NDDOAtom withNewParams(NDDOParams np);

	/**
	 * Returns a brand new NDDOAtom object. Everything is pass-by-value.
	 *
	 * @param coordinates New coordinates.
	 * @return New NDDOAtom.
	 */
	public abstract NDDOAtom withNewCoords(double[] coordinates);

	public NDDOOrbital[] getOrbitals() {
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

	public NDDOOrbital s() {
		return this.orbitals[0];
	}

	public double V(NDDOOrbital a, NDDOOrbital b) {
		return -this.atomProperties.getQ() * State.nom.getG(a, b, this.s(), this.s());
	}

	public double Vderiv(NDDOOrbital a, NDDOOrbital b, int tau) {
		return -this.atomProperties.getQ() * State.nom.getGderiv(a, b, this.s(), this.s(), tau);
	}

	public double Vderiv2(NDDOOrbital a, NDDOOrbital b, int tau1, int tau2) {
		return -this.atomProperties.getQ() * State.nom.getGderiv2(a, b, this.s(), this.s(), tau1, tau2);
	}

	public double VParamDeriv(NDDOOrbital a, NDDOOrbital b, int num, int type) {
		return -this.atomProperties.getQ() * State.nom.getGderiv(this.s(), this.s(), a, b, num, type);
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

	public double p1Deriv(int type) {
		double D1deriv = D1Deriv(type);

		if (getAtomProperties().getZ() == 1) {
			return 0;
		}

		return -p1 * p1 * D1 / (p1 * p1 * p1 -
				Math.pow(D1 * D1 + p1 * p1, 1.5)) * D1deriv;
	}

	public double p2Deriv(int type) {
		if (type == 0) {
			return 0;
		}

		double D2deriv = D2Deriv(type);

		double F1 = 2 * D2 * (Math.pow(D2 * D2 + p2 * p2, -1.5) -
				Math.pow(2 * D2 * D2 + p2 * p2, -1.5));

		double F2 = p2 * (2 * Math.pow(D2 * D2 + p2 * p2, -1.5) -
				Math.pow(2 * D2 * D2 + p2 * p2, -1.5)) -
				1 / (p2 * p2);

		return -F1 / F2 * D2deriv;

	}

	public double D1Deriv(int type) {
		double zetas = getParams().getZetas();

		double zetap = getParams().getZetap();

		double zeta = 0;

		switch (type) {

			case 0:
				zeta = zetap;
				break;
			case 1:
				zeta = zetas;
				break;
			default:
				zeta = 0;
		}

		return (2 * getAtomProperties().getPeriod() + 1) / Math.sqrt(3) *
				(4 * zeta * (0.5 + getAtomProperties().getPeriod()) *
						Math.pow(4 * zetas * zetap,
								-0.5 + getAtomProperties().getPeriod()) /
						Math.pow(zetas + zetap,
								2 + 2 * getAtomProperties().getPeriod())
						- (2 + 2 * getAtomProperties().getPeriod()) *
						Math.pow(4 * zetas * zetap,
								0.5 + getAtomProperties().getPeriod()) /
						Math.pow(zetas + zetap,
								3 + 2 * getAtomProperties().getPeriod()));


	}

	public double D2Deriv(int type) {
		if (type == 0) {
			return 0;
		}

		return -1 / getParams().getZetap() * D2;
	}

	public abstract double crf(NDDOAtom b);

	public abstract double crfDeriv(NDDOAtom b, int tau);

	public abstract double crfDeriv2(NDDOAtom b, int tau1, int tau2);

	public abstract double crfAlphaDeriv(NDDOAtom b, int num);

	@Override
	public abstract NDDOAtom clone();
}
