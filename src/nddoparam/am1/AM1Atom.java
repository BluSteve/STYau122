package nddoparam.am1;

import nddoparam.NDDO6G;
import nddoparam.NDDOAtom;
import nddoparam.NDDODerivative;
import scf.AtomProperties;
import scf.GTO;

import java.util.Objects;

public class AM1Atom extends NDDOAtom {
	private AM1Params mp;

	public AM1Atom(AtomProperties atomProperties, double[] coordinates, AM1Params mp) {
		super(atomProperties, coordinates, mp);
		this.mp = mp.clone();
	}

	// assign new coords
	public AM1Atom(AM1Atom a, double[] coords) {
		this(a.atomProperties, coords.clone(), a.mp.clone());
	}

	public AM1Atom(AM1Atom a) { // copy constructor
		this(a.atomProperties, a.coordinates.clone(), a.mp.clone());
	}

	// assign new params
	public AM1Atom(AM1Atom a, AM1Params mp) {
		this(a.atomProperties, a.coordinates.clone(), mp.clone());
	}


	private static double getf(AM1Atom a, AM1Atom b, double R) {
		return 1 + R / bohr * Math.exp(-a.mp.getAlpha() * R / bohr) +
				Math.exp(-b.mp.getAlpha() * R / bohr);
	}

	private static double getfPrime(AM1Atom a, AM1Atom b, double R, int tau) {
		return (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr) *
				Math.exp(-a.mp.getAlpha() * R / bohr)
				- a.mp.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) /
				(bohr * bohr) * Math.exp(-a.mp.getAlpha() * R / bohr)
				- b.mp.getAlpha() / bohr *
				(a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
				Math.exp(-b.mp.getAlpha() * R / bohr);
	}

	private static double F(AM1Atom a, double R) {
		return a.getParams().getK1() * Math.exp(
				-a.getParams().getL1() * (R - a.getParams().getM1()) *
						(R - a.getParams().getM1())) + a.getParams().getK2() * Math.exp(
				-a.getParams().getL2() * (R - a.getParams().getM2()) *
						(R - a.getParams().getM2())) + a.getParams().getK3() * Math.exp(
				-a.getParams().getL3() * (R - a.getParams().getM3()) *
						(R - a.getParams().getM3())) + a.getParams().getK4() * Math.exp(
				-a.getParams().getL4() * (R - a.getParams().getM4()) *
						(R - a.getParams().getM4()));

	}

	private static double Fderiv(AM1Atom a, double r, double val) {
		return gaussderiv(a.getParams().getK1(), a.getParams().getL1(),
				a.getParams().getM1(), r, val) +
				gaussderiv(a.getParams().getK2(), a.getParams().getL2(),
						a.getParams().getM2(), r, val) +
				gaussderiv(a.getParams().getK3(), a.getParams().getL3(),
						a.getParams().getM3(), r, val) +
				gaussderiv(a.getParams().getK4(), a.getParams().getL4(),
						a.getParams().getM4(), r, val);
	}

	private static double gaussderiv(double K, double L, double M, double R,
									 double val) {
		double r = R / bohr;
		double returnval = K * Math.exp(-L * (r - M) * (r - M));

		returnval = -returnval * 2 * L * (r - M) * val / (R * bohr);

		return returnval;
	}

	@Override
	public AM1Params getParams() { // TODO this is a workaround
		if (Objects.isNull(mp)) {
			return (AM1Params) super.np.clone();
		}
		return mp.clone();
	}

	@Override
	public double crf(NDDOAtom c) {
		AM1Atom b = (AM1Atom) c;
		double f;
		double R = GTO.R(this.getCoordinates(), b.getCoordinates());
		if ((this.atomProperties.getZ() == 7 || this.atomProperties.getZ() == 8) &&
				b.atomProperties.getZ() == 1) {
			f = getf(this, b, R);
		}
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) &&
				this.atomProperties.getZ() == 1) {
			f = getf(b, this, R);
		}
		else {
			f = 1 + Math.exp(-b.mp.getAlpha() * R / bohr) +
					Math.exp(-this.mp.getAlpha() * R / bohr);
		}

		return f * this.atomProperties.getQ() * b.atomProperties.getQ() *
				NDDO6G.getG(this.s(), this.s(), b.s(), b.s()) +
				(F(this, R) + F(b, R)) * this.getAtomProperties().getQ() *
						b.getAtomProperties().getQ() / R;
	}

	@Override
	public double crfDeriv(NDDOAtom c, int tau) {
		AM1Atom b = (AM1Atom) c;
		double R = GTO.R(this.getCoordinates(), b.getCoordinates());

		double r = R / bohr;
		double f;
		double fprime;
		if ((this.atomProperties.getZ() == 7 || this.atomProperties.getZ() == 8) &&
				b.atomProperties.getZ() == 1) {
			f = getf(this, b, R);
			fprime = getfPrime(this, b, R, tau);
		}
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) &&
				this.atomProperties.getZ() == 1) {
			f = getf(b, this, R);
			fprime = -getfPrime(b, this, R, tau);
		}
		else {
			f = 1 + Math.exp(-b.mp.getAlpha() * r) + Math.exp(-this.mp.getAlpha() * r);
			fprime = -b.mp.getAlpha() / bohr *
					(this.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
					Math.exp(-b.mp.getAlpha() * r)
					- this.mp.getAlpha() / bohr *
					(this.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
					Math.exp(-this.mp.getAlpha() * r);
		}

		double reciprocalderiv =
				-this.atomProperties.getQ() * b.atomProperties.getQ() / (r * r) *
						(this.getCoordinates()[tau] - b.getCoordinates()[tau]) /
						(R * bohr);

		double gaussiancontribution =
				(Fderiv(this, R, this.getCoordinates()[tau] - b.getCoordinates()[tau]) +
						Fderiv(b, R,
								this.getCoordinates()[tau] - b.getCoordinates()[tau])) *
						this.atomProperties.getQ() * b.atomProperties.getQ() / r +
						reciprocalderiv * (F(this, r) + F(b, r));

		double returnval =
				fprime * this.atomProperties.getQ() * b.atomProperties.getQ() *
				NDDO6G.getG(this.s(), this.s(), b.s(), b.s()) +
				f * this.atomProperties.getQ() * b.atomProperties.getQ() *
						NDDODerivative.getGderiv(this.s(), this.s(), b.s(), b.s(), tau) +
				gaussiancontribution;

		return returnval;
	}

	public double crfDeriv2(NDDOAtom c, int tau1, int tau2) {
		return 0;
	}

	@Override
	public double crfParamDeriv(NDDOAtom b, int num) {//todo ree
		return 0;
	}
}
