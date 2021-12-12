package nddo.am1;

import nddo.NDDO6G;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.geometry.GeometryDerivative;
import scf.AtomProperties;
import scf.GTO;

public class AM1Atom extends NDDOAtom {
	public static final int[] HParamNums =
			{0, 1, 3, 5, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18};
	public static final int[] NParamNums =
			{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18};
	public static final int[] OParamNums =
			{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 16, 17};
	public static final int[] CParamNums =
			{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18, 19};
	private final AM1Params mp;

	public AM1Atom(AtomProperties atomProperties, double[] coordinates, AM1Params mp) {
		super(atomProperties, coordinates, mp);
		this.mp = mp.clone();
	}

	@Override
	public NDDOAtom withNewParams(NDDOParams np) {
		return new AM1Atom(atomProperties, coordinates, (AM1Params) np); // todo bandaid cos of coming refactor
	}

	@Override
	public NDDOAtom withNewCoords(double[] coordinates) {
		return new AM1Atom(atomProperties, coordinates, mp);
	}

	@Override
	public NDDOAtom clone() {
		return new AM1Atom(atomProperties, coordinates, mp);
	}

	private static double getf(AM1Atom a, AM1Atom b, double R) {
		return 1 + R / bohr * Math.exp(-a.mp.getAlpha() * R / bohr) + Math.exp(-b.mp.getAlpha() * R / bohr);
	}

	private static double getfPrime(AM1Atom a, AM1Atom b, double R, int tau) {
		return (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr) *
				Math.exp(-a.mp.getAlpha() * R / bohr) -
				a.mp.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (bohr * bohr) *
						Math.exp(-a.mp.getAlpha() * R / bohr) -
				b.mp.getAlpha() / bohr * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
						Math.exp(-b.mp.getAlpha() * R / bohr);
	}

	private static double gaussderiv(double K, double L, double M, double R, double val) {
		double r = R / bohr;
		double returnval = K * Math.exp(-L * (r - M) * (r - M));

		returnval = -returnval * 2 * L * (r - M) * val / (R * bohr);

		return returnval;
	}

	private double Fderiv(double r, double val) {
		return gaussderiv(mp.getK1(), mp.getL1(), mp.getM1(), r, val) +
				gaussderiv(mp.getK2(), mp.getL2(), mp.getM2(), r, val) +
				gaussderiv(mp.getK3(), mp.getL3(), mp.getM3(), r, val) +
				gaussderiv(mp.getK4(), mp.getL4(), mp.getM4(), r, val);
	}

	private double F(double R) {
		return mp.getK1() * Math.exp(-mp.getL1() * (R - mp.getM1()) * (R - mp.getM1())) +
				mp.getK2() * Math.exp(-mp.getL2() * (R - mp.getM2()) * (R - mp.getM2())) +
				mp.getK3() * Math.exp(-mp.getL3() * (R - mp.getM3()) * (R - mp.getM3())) +
				mp.getK4() * Math.exp(-mp.getL4() * (R - mp.getM4()) * (R - mp.getM4()));
	}

	@Override
	public double crf(NDDOAtom c) {
		AM1Atom b = (AM1Atom) c;
		double f;
		double R = GTO.R(this.getCoordinates(), b.getCoordinates());

		if ((this.atomProperties.getZ() == 7 || this.atomProperties.getZ() == 8) && b.atomProperties.getZ() == 1)
			f = getf(this, b, R);
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && this.atomProperties.getZ() == 1)
			f = getf(b, this, R);
		else f = 1 + Math.exp(-b.mp.getAlpha() * R / bohr) + Math.exp(-this.mp.getAlpha() * R / bohr);

		return f * this.atomProperties.getQ() * b.atomProperties.getQ() *
				NDDO6G.getG(this.s(), this.s(), b.s(), b.s()) +
				(this.F(R) + b.F(R)) * this.getAtomProperties().getQ() * b.getAtomProperties().getQ() / R;
	}

	@Override
	public double crfDeriv(NDDOAtom c, int tau) {
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
			f = 1 + Math.exp(-b.mp.getAlpha() * r) + Math.exp(-this.mp.getAlpha() * r);
			fprime = -b.mp.getAlpha() / bohr * (this.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
					Math.exp(-b.mp.getAlpha() * r) -
					this.mp.getAlpha() / bohr * (this.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
							Math.exp(-this.mp.getAlpha() * r);
		}

		double reciprocalderiv = -this.atomProperties.getQ() * b.atomProperties.getQ() / (r * r) *
				(this.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr);

		double gaussiancontribution = (this.Fderiv(R, this.getCoordinates()[tau] - b.getCoordinates()[tau]) +
				b.Fderiv(R, this.getCoordinates()[tau] - b.getCoordinates()[tau])) * this.atomProperties.getQ() *
				b.atomProperties.getQ() / r + reciprocalderiv * (this.F(r) + b.F(r));

		return fprime * this.atomProperties.getQ() * b.atomProperties.getQ() *
				NDDO6G.getG(this.s(), this.s(), b.s(), b.s()) +
				f * this.atomProperties.getQ() * b.atomProperties.getQ() *
						GeometryDerivative.getGderiv(this.s(), this.s(), b.s(), b.s(), tau) + gaussiancontribution;
	}

	@Override
	public double crfDeriv2(NDDOAtom c, int tau1, int tau2) {
		return 0;
	}

	@Override
	public double crfAlphaDeriv(NDDOAtom b, int num) {//todo ree
		return 0;
	}
}
