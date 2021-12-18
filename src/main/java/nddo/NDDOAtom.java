package nddo;

import nddo.structs.AtomProperties;

public abstract class NDDOAtom { // todo make this interface
	public final double p0, p1, p2, D1, D2;
	protected final double[] coordinates;
	protected final AtomProperties atomProperties;
	protected final NDDOParams np;
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

		if (this.atomProperties.getZ() == 1) {
			this.p1 = 0;
			this.p2 = 0;
			this.D1 = 0;
			this.D2 = 0;
		}
		else {
			this.D1 = D1();
			this.D2 = D2();
			this.p1 = p1();
			this.p2 = p2();
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
	 * Returns a brand new NDDOAtom object. Deep copied!
	 *
	 * @param np New NDDOParams.
	 * @return New NDDOAtom.
	 */
	public abstract NDDOAtom withNewParams(NDDOParams np);

	/**
	 * Returns a brand new NDDOAtom object. Not cloned!
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

	protected NDDOOrbital s() {
		return this.orbitals[0];
	}

	public double V(NDDOOrbital a, NDDOOrbital b) {
		return -this.atomProperties.getQ() * State.nom.G(a, b, this.s(), this.s());
	}

	public double Vgd(NDDOOrbital a, NDDOOrbital b, int tau) {
		return -this.atomProperties.getQ() * State.nom.Ggd(a, b, this.s(), this.s(), tau);
	}

	public double Vg2d(NDDOOrbital a, NDDOOrbital b, int tau1, int tau2) {
		return -this.atomProperties.getQ() * State.nom.Gg2d(a, b, this.s(), this.s(), tau1, tau2);
	}

	public double Vpd(NDDOOrbital a, NDDOOrbital b, int num, int type) {
		return -this.atomProperties.getQ() * State.nom.Gpd(this.s(), this.s(), a, b, num, type);
	}

	protected double p0() {
		return Constants.eV / (2 * np.getGss());
	}

	protected double p1() {
		double guess = 0;
		double newguess = 0.5 * Math.pow(D1 * D1 * Constants.eV / np.getHsp(), 1.0 / 3);
		while (Math.abs(guess - newguess) > 1E-12) {
			guess = newguess;
			double f = 1 / guess - 1 / Math.sqrt(guess * guess + D1 * D1) - 4 * np.getHsp() / Constants.eV;
			double fprime = -1 / (guess * guess) + guess / Math.pow(guess * guess + D1 * D1, 1.5);
			newguess = guess - f / fprime;
		}
		return newguess;
	}

	public double p1pd(int type) {
		double D1deriv = D1pd(type);

		if (getAtomProperties().getZ() == 1) {
			return 0;
		}

		return -p1 * p1 * D1 / (p1 * p1 * p1 -
				Math.pow(D1 * D1 + p1 * p1, 1.5)) * D1deriv;
	}

	public double p1p2d(int type) {
		double D1 = this.D1;
		double p1 = this.p1;

		double D1sderiv = D1pd(0);
		double D1pderiv = p1pd(1);
		double p1sderiv = p1pd(0);
		double p1pderiv = p1pd(1);

		double D1ssderiv2 = D1p2d(0);
		double D1spderiv2 = D1p2d(1);
		double D1ppderiv2 = D1p2d(2);

		double num0 = Math.sqrt(D1 * D1 + p1 * p1);
		double num1 = Math.pow(p1 * p1 + D1 * D1, 1.5);

		switch (type) {
			case 0:
				return -(3 * (p1 * p1 - p1 * num0) * p1sderiv * p1sderiv +
						D1 * (2 * p1 - 3 * num0) * D1sderiv * p1sderiv + p1 * p1 * D1sderiv * D1sderiv +
						p1 * p1 * D1 * D1ssderiv2) / (p1 * p1 * p1 - num1);
			case 1:
				return -(3 * (p1 * p1 - p1 * num0) * p1sderiv * p1pderiv + D1 * 2 * p1 * D1sderiv * p1pderiv -
						3 * D1 * num0 * D1pderiv * p1sderiv + p1 * p1 * D1sderiv * D1pderiv +
						p1 * p1 * D1 * D1spderiv2) / (p1 * p1 * p1 - num1);
			case 2:
				return -(3 * (p1 * p1 - p1 * num0) * p1pderiv * p1pderiv +
						D1 * (2 * p1 - 3 * num0) * D1pderiv * p1pderiv + p1 * p1 * D1pderiv * D1pderiv +
						p1 * p1 * D1 * D1ppderiv2) / (p1 * p1 * p1 - num1);

		}

		return 0;
	}

	protected double p2() {
		double guess = 0;
		double newguess = 0.5;
		while (Math.abs(guess - newguess) > 1E-12) {
			guess = newguess;
			double f = 1 / guess + 1 / Math.sqrt(guess * guess + 2 * D2 * D2) -
					2 / Math.sqrt(guess * guess + D2 * D2) - 4 * (np.getGpp() - np.getGp2()) / Constants.eV;
			double fprime = -1 / (guess * guess) - guess / Math.pow(guess * guess + 2 * D2 * D2, 1.5) +
					2 * guess / Math.pow(guess * guess + D2 * D2, 1.5);
			newguess = guess - f / fprime;
		}
		return newguess;
	}

	public double p2pd(int type) {
		if (type == 0) {
			return 0;
		}

		double D2deriv = D2pd(type);

		double F1 = 2 * D2 * (Math.pow(D2 * D2 + p2 * p2, -1.5) -
				Math.pow(2 * D2 * D2 + p2 * p2, -1.5));

		double F2 = p2 * (2 * Math.pow(D2 * D2 + p2 * p2, -1.5) -
				Math.pow(2 * D2 * D2 + p2 * p2, -1.5)) -
				1 / (p2 * p2);

		return -F1 / F2 * D2deriv;

	}

	public double p2p2d(int type) {
		if (type != 2) {
			return 0;
		}

		double D2 = this.D2;
		double p2 = this.p2;
		double D2pderiv = D2pd(1);
		double p2pderiv = p2pd(1);
		double D2ppderiv2 = D2p2d(2);

		double num0 = D2 * D2 * 2 + p2 * p2;
		double num1 = Math.pow(D2 * D2 + p2 * p2, -1.5);
		double num2 = Math.pow(num0, -1.5);
		double num3 = Math.pow(D2 * D2 + p2 * p2, -2.5);
		double num4 = (2 * D2 * D2pderiv + p2 * p2pderiv) * Math.pow(num0, -2.5);

		double F1 = 2 * D2 * (num1 - num2);
		double F2 = p2 * (2 * num1 - num2) - 1 / (p2 * p2);

		double F1deriv = 2 * D2pderiv * (num1 - num2) - 6 * D2 * ((D2 * D2pderiv + p2 * p2pderiv) * num3 - num4);
		double F2deriv = p2pderiv * (2 * num1 - num2) + 2 / (p2 * p2 * p2) * p2pderiv -
				3 * p2 * (2 * (D2 * D2pderiv + p2 * p2pderiv) * num3 - num4);

		return -F1 / F2 * D2ppderiv2 - F1deriv / F2 * D2pderiv + F1 / (F2 * F2) * F2deriv * D2pderiv;
	}

	protected double D1() {
		return (2 * atomProperties.getPeriod() + 1) / Math.sqrt(3) *
				Math.pow(4 * np.getZetas() * np.getZetap(), atomProperties.getPeriod() + 0.5) /
				Math.pow(np.getZetas() + np.getZetap(), 2 * atomProperties.getPeriod() + 2);
	}

	public double D1pd(int type) {
		double zetas = getParams().getZetas();

		double zetap = getParams().getZetap();

		double zeta;

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

	public double D1p2d(int type) {
		double zetas = getParams().getZetas();
		double zetap = getParams().getZetap();
		int n = getAtomProperties().getPeriod();

		double num0 = Math.pow(4 * zetas * zetap, -1.5 + n);
		double num1 = Math.pow(zetas + zetap, 2 * n + 2);
		double num2 = Math.pow(4 * zetas * zetap, -0.5 + n);
		double num3 = Math.pow(zetas + zetap, 2 * n + 3);
		double num4 = (2 * n + 2) * (2 * n + 3) * Math.pow(4 * zetas * zetap, 0.5 + n)
				/ Math.pow(zetas + zetap, 4 + 2 * n);

		switch (type) {
			case 0:
				return (2 * getAtomProperties().getPeriod() + 1) / Math.sqrt(3) *
						(16 * zetap * zetap * (n + 0.5) * (n - 0.5) * num0 / num1 -
								8 * zetap * (n + 0.5) * (2 * n + 2) * num2 / num3 + num4);
			case 1:
				return (2 * getAtomProperties().getPeriod() + 1) / Math.sqrt(3) *
						(4 * (n + 0.5) * Math.pow(4 * zetas * zetap, n - 0.5) / num1 +
								16 * zetas * zetap * (n + 0.5) * (n - 0.5) * num0 / num1 -
								4 * (zetas + zetap) * (n + 0.5) * (2 * n + 2) * num2 / num3 + num4);
			case 2:
				return (2 * getAtomProperties().getPeriod() + 1) / Math.sqrt(3) *
						(16 * zetas * zetas * (n + 0.5) * (n - 0.5) * num0 / num1 -
								8 * zetas * (n + 0.5) * (2 * n + 2) * num2 / num3 + num4);

		}

		return 0;
	}

	protected double D2() {
		return 1 / np.getZetap() *
				Math.sqrt((2 * atomProperties.getPeriod() + 1) * (2 * atomProperties.getPeriod() + 2) / 20.0);
	}

	public double D2pd(int type) {
		if (type == 0) {
			return 0;
		}

		return -1 / getParams().getZetap() * D2;
	}

	public double D2p2d(int type) {
		if (type == 2) {
			return 2 * this.D2 / (getParams().getZetap() * getParams().getZetap());
		}

		return 0;
	}

	public abstract double crf(NDDOAtom b);

	public abstract double crfDeriv(NDDOAtom b, int tau);

	public abstract double crfDeriv2(NDDOAtom b, int tau1, int tau2);

	public abstract double crfAlphaDeriv(NDDOAtom b, int num);

	@Override
	public abstract NDDOAtom clone();
}
