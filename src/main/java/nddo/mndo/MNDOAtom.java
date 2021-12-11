package nddo.mndo;

import nddo.NDDO6G;
import nddo.NDDOAtom;
import nddo.geometry.GeometryDerivative;
import scf.AtomHandler;
import scf.AtomProperties;
import scf.GTO;

import java.io.IOException;

public class MNDOAtom extends NDDOAtom {
	private final MNDOParams mp;

	public MNDOAtom(AtomProperties atomProperties, double[] coordinates, MNDOParams mp) {
		super(atomProperties, coordinates, mp);
		this.mp = mp.clone();
	}

	private static double getf(MNDOAtom a, MNDOAtom b, double R) {
		return 1 + R / bohr * Math.exp(-a.mp.getAlpha() * R / bohr) + Math.exp(-b.mp.getAlpha() * R / bohr);
	}

	private static double getfPrime(MNDOAtom a, MNDOAtom b, double R, int tau) {
		return (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr) *
				Math.exp(-a.mp.getAlpha() * R / bohr) -
				a.mp.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (bohr * bohr) *
						Math.exp(-a.mp.getAlpha() * R / bohr) -
				b.mp.getAlpha() / bohr * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
						Math.exp(-b.mp.getAlpha() * R / bohr);
	}

	@Override
	public double crf(NDDOAtom c) {
		MNDOAtom b = (MNDOAtom) c;
		double f;
		double R = GTO.R(this.getCoordinates(), b.getCoordinates());

		if ((this.atomProperties.getZ() == 7 || this.atomProperties.getZ() == 8) && b.atomProperties.getZ() == 1)
			f = getf(this, b, R);
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && this.atomProperties.getZ() == 1)
			f = getf(b, this, R);
		else f = 1 + Math.exp(-b.mp.getAlpha() * R / bohr) + Math.exp(-this.mp.getAlpha() * R / bohr);

		return f * this.atomProperties.getQ() * b.atomProperties.getQ() * NDDO6G.getG(this.s(), this.s(), b.s(),
				b.s());
	}

	@Override
	public double crfDeriv(NDDOAtom c, int tau) {
		MNDOAtom b = (MNDOAtom) c;
		double f;
		double fprime;
		double R = GTO.R(this.getCoordinates(), b.getCoordinates());
		double r = R / bohr;

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

		return fprime * this.atomProperties.getQ() * b.atomProperties.getQ() *
				NDDO6G.getG(this.s(), this.s(), b.s(), b.s()) +
				f * this.atomProperties.getQ() * b.atomProperties.getQ() *
						GeometryDerivative.getGderiv(this.s(), this.s(), b.s(), b.s(), tau);
	}

	@Override
	public double crfDeriv2(NDDOAtom c, int tau1, int tau2) {/*TODO*/
		double orig = this.crfDeriv(c, tau1);
		double[] coords = this.getCoordinates().clone();

		coords[tau2] = coords[tau2] + 1E-7;
		NDDOAtom perturbed = new MNDOAtom(this.atomProperties, coords, this.mp);

		double newval = perturbed.crfDeriv(c, tau1);

		return 1E7 * (newval - orig);
	}

	@Override
	public double crfParamDeriv(NDDOAtom b, int num) {
		MNDOAtom c = (MNDOAtom) b;
		double R = GTO.R(this.getCoordinates(), b.getCoordinates());
		double val = this.atomProperties.getQ() * c.atomProperties.getQ() *
				NDDO6G.getG(this.s(), this.s(), b.s(), b.s());

		double returnval = 0;

		if (num == 0 || num == 2) returnval += val * -R / bohr * Math.exp(-this.mp.getAlpha() * R / bohr);
		if (num == 1 || num == 2) returnval += val * -R / bohr * Math.exp(-c.mp.getAlpha() * R / bohr);

		return returnval;
	}
}
