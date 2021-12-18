package examples.am1;

import nddo.NDDOParams;
import nddo.State;
import nddo.defaults.NDDO6G;
import nddo.defaults.NDDOAtomBasic;
import nddo.scf.GTO;
import nddo.structs.AtomProperties;
import nddo.structs.OrbitalProperties;

import static nddo.Constants.bohr;

public class AM1Atom extends NDDOAtomBasic {
	public AM1Atom(AtomProperties atomProperties, double[] coordinates, NDDOParams np) {
		super(atomProperties, coordinates, np);

		OrbitalProperties[] orbitalProperties = this.atomProperties.getOrbitals();
		orbitals = new NDDO6G[orbitalProperties.length];

		for (int x = 0; x < orbitals.length; x++) {
			switch (orbitalProperties[x].type) {
				case "s":
					orbitals[x] = new NDDO6G(this, orbitalProperties[x], np.getZetas(), np.getBetas(), np.getUss());
					break;
				case "p":
					orbitals[x] = new NDDO6G(this, orbitalProperties[x], np.getZetap(), np.getBetap(), np.getUpp());
					break;
			}
		}
	}

	private static double getf(AM1Atom a, AM1Atom b, double R) {
		return 1 + R / bohr * Math.exp(-a.np.getAlpha() * R / bohr) + Math.exp(-b.np.getAlpha() * R / bohr);
	}

	private static double getfPrime(AM1Atom a, AM1Atom b, double R, int tau) {
		return (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr) *
				Math.exp(-a.np.getAlpha() * R / bohr) -
				a.np.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (bohr * bohr) *
						Math.exp(-a.np.getAlpha() * R / bohr) -
				b.np.getAlpha() / bohr * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
						Math.exp(-b.np.getAlpha() * R / bohr);
	}

	private static double gaussderiv(double K, double L, double M, double R, double val) {
		double r = R / bohr;
		double returnval = K * Math.exp(-L * (r - M) * (r - M));

		returnval = -returnval * 2 * L * (r - M) * val / (R * bohr);

		return returnval;
	}

	@Override
	public NDDOAtomBasic withNewParams(NDDOParams np) {
		return new AM1Atom(atomProperties, coordinates, np); // todo bandaid cos of coming refactor
	}

	@Override
	public NDDOAtomBasic withNewCoords(double[] coordinates) {
		return new AM1Atom(atomProperties, coordinates, np);
	}

	private double Fderiv(double r, double val) {
		return gaussderiv(np.aget(0), np.aget(4), np.aget(8), r, val) +
				gaussderiv(np.aget(1), np.aget(5), np.aget(9), r, val) +
				gaussderiv(np.aget(2), np.aget(6), np.aget(10), r, val) +
				gaussderiv(np.aget(3), np.aget(7), np.aget(11), r, val);
	}

	private double F(double R) {
		return np.aget(0) * Math.exp(-np.aget(4) * (R - np.aget(8)) * (R - np.aget(8))) +
				np.aget(1) * Math.exp(-np.aget(5) * (R - np.aget(9)) * (R - np.aget(9))) +
				np.aget(2) * Math.exp(-np.aget(6) * (R - np.aget(10)) * (R - np.aget(10))) +
				np.aget(3) * Math.exp(-np.aget(7) * (R - np.aget(11)) * (R - np.aget(11)));
	}

	@Override
	public double crf(NDDOAtomBasic c) {
		AM1Atom b = (AM1Atom) c;
		double f;
		double R = GTO.R(this.getCoordinates(), b.getCoordinates());

		if ((this.atomProperties.getZ() == 7 || this.atomProperties.getZ() == 8) && b.atomProperties.getZ() == 1)
			f = getf(this, b, R);
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && this.atomProperties.getZ() == 1)
			f = getf(b, this, R);
		else f = 1 + Math.exp(-b.np.getAlpha() * R / bohr) + Math.exp(-this.np.getAlpha() * R / bohr);

		return f * this.atomProperties.getQ() * b.atomProperties.getQ() *
				State.nom.G(this.s(), this.s(), b.s(), b.s()) +
				(this.F(R) + b.F(R)) * this.getAtomProperties().getQ() * b.getAtomProperties().getQ() / R;
	}

	@Override
	public double crfgd(NDDOAtomBasic c, int tau) {
		AM1Atom b = (AM1Atom) c;
		double R = GTO.R(this.getCoordinates(), b.getCoordinates());
		double r = R / bohr;
		double f;
		double fprime;

		if ((this.atomProperties.getZ() == 7 || this.atomProperties.getZ() == 8) && b.atomProperties.getZ() == 1) {
			f = getf(this, b, R);
			fprime = getfPrime(this, b, R, tau);
		}
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && this.atomProperties.getZ() == 1) {
			f = getf(b, this, R);
			fprime = -getfPrime(b, this, R, tau);
		}
		else {
			f = 1 + Math.exp(-b.np.getAlpha() * r) + Math.exp(-this.np.getAlpha() * r);
			fprime = -b.np.getAlpha() / bohr * (this.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
					Math.exp(-b.np.getAlpha() * r) -
					this.np.getAlpha() / bohr * (this.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
							Math.exp(-this.np.getAlpha() * r);
		}

		double reciprocalderiv = -this.atomProperties.getQ() * b.atomProperties.getQ() / (r * r) *
				(this.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr);

		double gaussiancontribution = (this.Fderiv(R, this.getCoordinates()[tau] - b.getCoordinates()[tau]) +
				b.Fderiv(R, this.getCoordinates()[tau] - b.getCoordinates()[tau])) * this.atomProperties.getQ() *
				b.atomProperties.getQ() / r + reciprocalderiv * (this.F(r) + b.F(r));

		return fprime * this.atomProperties.getQ() * b.atomProperties.getQ() *
				State.nom.G(this.s(), this.s(), b.s(), b.s()) +
				f * this.atomProperties.getQ() * b.atomProperties.getQ() *
						State.nom.Ggd(this.s(), this.s(), b.s(), b.s(), tau) + gaussiancontribution;
	}

	@Override
	public double crfg2d(NDDOAtomBasic c, int tau1, int tau2) {
		return 0;
	}

	@Override
	public double crfpd(NDDOAtomBasic b, int num) {//todo ree
		return 0;
	}

	@Override
	public NDDOAtomBasic copy() {
		return new AM1Atom(atomProperties, coordinates, np);
	}
}
