package nddo.defaults;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.State;
import nddo.structs.AtomProperties;
import tools.Pow;

// some methods that might come in handy if you're doing your own NDDO impl.
public abstract class NDDOAtomBasic<T extends NDDOAtom> implements NDDOAtom<T, NDDO6G> {
	protected final double p0, p1, p2, D1, D2;
	protected final double[] coordinates;
	protected final AtomProperties atomProperties;
	protected final NDDOParams np;
	protected NDDO6G[] orbitals;

	public NDDOAtomBasic(AtomProperties atomProperties, double[] coordinates, NDDOParams np) {
		this.atomProperties = atomProperties;
		this.coordinates = coordinates;
		this.np = np;

		this.p0 = findp0();

		if (this.atomProperties.getZ() == 1) {
			this.p1 = 0;
			this.p2 = 0;
			this.D1 = 0;
			this.D2 = 0;
		}
		else {
			this.D1 = findD1();
			this.D2 = findD2();
			this.p1 = findp1();
			this.p2 = findp2();
		}
	}

	@Override
	public double p0() {
		return p0;
	}

	@Override
	public double p1() {
		return p1;
	}

	@Override
	public double p2() {
		return p2;
	}

	@Override
	public double D1() {
		return D1;
	}

	@Override
	public double D2() {
		return D2;
	}

	/**
	 * Get coordinates object with no cloning. Can be modified in place!
	 *
	 * @return Original coordinates array object.
	 */
	@Override
	public double[] getCoordinates() {
		return this.coordinates;
	}

	@Override
	public AtomProperties getAtomProperties() {
		return this.atomProperties;
	}

	/**
	 * Returns a brand new NDDOAtom object. Not copied!
	 *
	 * @param np New NDDOParams.
	 * @return New NDDOAtom.
	 */
	@Override
	public abstract T withNewParams(NDDOParams np);

	/**
	 * Returns a brand new NDDOAtom object. Not copied!
	 *
	 * @param coordinates New coordinates.
	 * @return New NDDOAtom.
	 */
	@Override
	public abstract T withNewCoords(double[] coordinates);

	@Override
	public NDDO6G[] getOrbitals() {
		return this.orbitals;
	}

	@Override
	public NDDOParams getParams() {
		return np;
	}

	@Override
	public NDDO6G s() {
		return this.orbitals[0];
	}

	@Override
	public double V(NDDO6G a, NDDO6G b) {
		return -this.atomProperties.getQ() * State.nom.G(a, b, this.s(), this.s());
	}

	@Override
	public double Vgd(NDDO6G a, NDDO6G b, int tau) {
		return -this.atomProperties.getQ() * State.nom.Ggd(a, b, this.s(), this.s(), tau);
	}

	@Override
	public double Vg2d(NDDO6G a, NDDO6G b, int tau1, int tau2) {
		return -this.atomProperties.getQ() * State.nom.Gg2d(a, b, this.s(), this.s(), tau1, tau2);
	}

	@Override
	public double Vpd(NDDO6G a, NDDO6G b, int num, int type) {
		return -this.atomProperties.getQ() * State.nom.Gpd(this.s(), this.s(), a, b, num, type);
	}

	@Override
	public double Vp2d(NDDO6G a, NDDO6G b, int num1, int type1, int num2, int type2) {
		return -this.atomProperties.getQ() * State.nom.Gp2d(this.s(), this.s(), a, b, num1, type1, num2, type2);
	}

	@Override
	public double Vpgd(NDDO6G a, NDDO6G b, int num, int type, int tau) {
		return -this.atomProperties.getQ() * State.nom.Gpgd(a, b, this.s(), this.s(), num, type, tau);

	}

	@Override
	public double Vp2gd(NDDO6G a, NDDO6G b, int num1, int type1, int num2, int type2, int tau) {
		return -this.atomProperties.getQ() * State.nom.Gp2gd(a, b, this.s(), this.s(), num1, type1, num2, type2, tau);

	}

	protected double findp0() {
		return Constants.eV / (2 * np.getGss());
	}

	protected double findp1() {
		double guess = 0;
		double newguess = 0.5 * Math.pow(D1 * D1 * Constants.eV / np.getHsp(), 1.0 / 3);
		while (Math.abs(guess - newguess) > 1E-12) {
			guess = newguess;
			double f = 1 / guess - 1 / Math.sqrt(guess * guess + D1 * D1) - 4 * np.getHsp() / Constants.eV;
			double fprime = -1 / (guess * guess) + guess / Pow.pow(guess * guess + D1 * D1, 1.5);
			newguess = guess - f / fprime;
		}
		return newguess;
	}

	@Override
	public double p1pd(int type) {
		double D1deriv = D1pd(type);

		if (atomProperties.getZ() == 1) {
			return 0;
		}

		return -p1 * p1 * D1 / (p1 * p1 * p1 -
				Pow.pow(D1 * D1 + p1 * p1, 1.5)) * D1deriv;
	}

	@Override
	public double p1p2d(int type) {
		double D1 = this.D1;
		double p1 = this.p1;

		double D1sderiv = D1pd(0);
		double D1pderiv = D1pd(1);
		double p1sderiv = p1pd(0);
		double p1pderiv = p1pd(1);
		double D1ssderiv2 = D1p2d(0);
		double D1spderiv2 = D1p2d(1);
		double D1ppderiv2 = D1p2d(2);

		double num0 = Math.sqrt(D1 * D1 + p1 * p1);
		double num1 = Pow.pow(p1 * p1 + D1 * D1, 1.5);
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

	protected double findp2() {
		double guess = 0;
		double newguess = 0.5;
		while (Math.abs(guess - newguess) > 1E-12) {
			guess = newguess;
			double a = guess * guess + 2 * D2 * D2;
			double f = 1 / guess + 1 / Math.sqrt(a) - 2 / Math.sqrt(guess * guess + D2 * D2) -
					4 * (np.getGpp() - np.getGp2()) / Constants.eV;
			double fprime = -1 / (guess * guess) - guess / Pow.pow(a, 1.5) +
					2 * guess / Pow.pow(guess * guess + D2 * D2, 1.5);
			newguess = guess - f / fprime;
		}
		return newguess;
	}

	@Override
	public double p2pd(int type) {
		if (type == 0) {
			return 0;
		}

		double D2deriv = D2pd(type);

		double pow = Pow.pow(D2 * D2 + p2 * p2, -1.5);
		double pow1 = Pow.pow(2 * D2 * D2 + p2 * p2, -1.5);

		double F1 = 2 * D2 * (pow - pow1);
		double F2 = p2 * (2 * pow - pow1) - 1 / (p2 * p2);

		return -F1 / F2 * D2deriv;
	}

	@Override
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
		double num1 = Pow.pow(D2 * D2 + p2 * p2, -1.5);
		double num2 = Pow.pow(num0, -1.5);
		double num3 = Pow.pow(D2 * D2 + p2 * p2, -2.5);
		double num4 = (2 * D2 * D2pderiv + p2 * p2pderiv) * Pow.pow(num0, -2.5);

		double F1 = 2 * D2 * (num1 - num2);
		double F2 = p2 * (2 * num1 - num2) - 1 / (p2 * p2);

		double F1deriv = 2 * D2pderiv * (num1 - num2) - 6 * D2 * ((D2 * D2pderiv + p2 * p2pderiv) * num3 - num4);
		double F2deriv = p2pderiv * (2 * num1 - num2) + 2 / (p2 * p2 * p2) * p2pderiv -
				3 * p2 * (2 * (D2 * D2pderiv + p2 * p2pderiv) * num3 - num4);

		return -F1 / F2 * D2ppderiv2 - F1deriv / F2 * D2pderiv + F1 / (F2 * F2) * F2deriv * D2pderiv;
	}

	protected double findD1() {
		return (2 * atomProperties.getPeriod() + 1) / Math.sqrt(3) *
				Pow.pow(4 * np.getZetas() * np.getZetap(), atomProperties.getPeriod() + 0.5) /
				Pow.pow(np.getZetas() + np.getZetap(), 2 * atomProperties.getPeriod() + 2);
	}

	@Override
	public double D1pd(int type) {
		double zetas = np.getZetas();

		double zetap = np.getZetap();

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

		return (2 * atomProperties.getPeriod() + 1) / Math.sqrt(3) *
				(4 * zeta * (0.5 + atomProperties.getPeriod()) *
						Pow.pow(4 * zetas * zetap, -0.5 + atomProperties.getPeriod()) /
						Pow.pow(zetas + zetap, 2 + 2 * atomProperties.getPeriod())
						- (2 + 2 * atomProperties.getPeriod()) *
						Pow.pow(4 * zetas * zetap, 0.5 + atomProperties.getPeriod()) /
						Pow.pow(zetas + zetap, 3 + 2 * atomProperties.getPeriod()));


	}

	@Override
	public double D1p2d(int type) {
		double zetas = np.getZetas();
		double zetap = np.getZetap();
		int n = atomProperties.getPeriod();

		double num0 = Pow.pow(4 * zetas * zetap, -1.5 + n);
		double num1 = Pow.pow(zetas + zetap, 2 * n + 2);
		double num2 = Pow.pow(4 * zetas * zetap, -0.5 + n);
		double num3 = Pow.pow(zetas + zetap, 2 * n + 3);
		double num4 = (2 * n + 2) * (2 * n + 3) * Pow.pow(4 * zetas * zetap, 0.5 + n)
				/ Pow.pow(zetas + zetap, 4 + 2 * n);

		switch (type) {
			case 0:
				return (2 * atomProperties.getPeriod() + 1) / Math.sqrt(3) *
						(16 * zetap * zetap * (n + 0.5) * (n - 0.5) * num0 / num1 -
								8 * zetap * (n + 0.5) * (2 * n + 2) * num2 / num3 + num4);
			case 1:
				return (2 * atomProperties.getPeriod() + 1) / Math.sqrt(3) *
						(4 * (n + 0.5) * Pow.pow(4 * zetas * zetap, n - 0.5) / num1 +
								16 * zetas * zetap * (n + 0.5) * (n - 0.5) * num0 / num1 -
								4 * (zetas + zetap) * (n + 0.5) * (2 * n + 2) * num2 / num3 + num4);
			case 2:
				return (2 * atomProperties.getPeriod() + 1) / Math.sqrt(3) *
						(16 * zetas * zetas * (n + 0.5) * (n - 0.5) * num0 / num1 -
								8 * zetas * (n + 0.5) * (2 * n + 2) * num2 / num3 + num4);
		}

		return 0;
	}

	protected double findD2() {
		return 1 / np.getZetap() *
				Math.sqrt((2 * atomProperties.getPeriod() + 1) * (2 * atomProperties.getPeriod() + 2) / 20.0);
	}

	@Override
	public double D2pd(int type) {
		if (type == 0) {
			return 0;
		}

		return -1 / np.getZetap() * D2;
	}

	@Override
	public double D2p2d(int type) {
		if (type == 2) {
			return 2 * this.D2 / (np.getZetap() * np.getZetap());
		}

		return 0;
	}
}
