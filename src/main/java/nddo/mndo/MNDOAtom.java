package nddo.mndo;

import nddo.NDDO6G;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.geometry.GeometryDerivative;
import nddo.structs.AtomProperties;
import nddo.scf.GTO;
import nddo.structs.OrbitalProperties;

import static nddo.Constants.bohr;
import static nddo.State.nom;

public class MNDOAtom extends NDDOAtom {
	public static final int[] T1ParamNums = {0, 1, 3, 5, 7};
	public static final int[] T2ParamNums = {0, 1, 2, 3, 4, 5, 6, 7};

	public MNDOAtom(AtomProperties atomProperties, double[] coordinates, NDDOParams np) {
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

	@Override
	public NDDOAtom withNewParams(NDDOParams np) {
		return new MNDOAtom(atomProperties, coordinates, np);
	}

	@Override
	public NDDOAtom withNewCoords(double[] coordinates) {
		return new MNDOAtom(atomProperties, coordinates, np);
	}

	@Override
	public NDDOAtom clone() {
		return new MNDOAtom(atomProperties, coordinates, np);
	}

	private static double getf(MNDOAtom a, MNDOAtom b, double R) {
		return 1 + R / bohr * Math.exp(-a.np.getAlpha() * R / bohr) + Math.exp(-b.np.getAlpha() * R / bohr);
	}

	private static double getfPrime(MNDOAtom a, MNDOAtom b, double R, int tau) {
		return (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr) *
				Math.exp(-a.np.getAlpha() * R / bohr) -
				a.np.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (bohr * bohr) *
						Math.exp(-a.np.getAlpha() * R / bohr) -
				b.np.getAlpha() / bohr * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
						Math.exp(-b.np.getAlpha() * R / bohr);
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
		else f = 1 + Math.exp(-b.np.getAlpha() * R / bohr) + Math.exp(-this.np.getAlpha() * R / bohr);

		return f * this.atomProperties.getQ() * b.atomProperties.getQ() * nom.getG(this.s(), this.s(), b.s(),
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
			f = 1 + Math.exp(-b.np.getAlpha() * r) + Math.exp(-this.np.getAlpha() * r);
			fprime = -b.np.getAlpha() / bohr * (this.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
					Math.exp(-b.np.getAlpha() * r) -
					this.np.getAlpha() / bohr * (this.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
							Math.exp(-this.np.getAlpha() * r);
		}

		return fprime * this.atomProperties.getQ() * b.atomProperties.getQ() *
				nom.getG(this.s(), this.s(), b.s(), b.s()) +
				f * this.atomProperties.getQ() * b.atomProperties.getQ() *
						nom.getGderiv(this.s(), this.s(), b.s(), b.s(), tau);
	}

	@Override
	public double crfDeriv2(NDDOAtom c, int tau1, int tau2) {/*TODO*/
		double orig = this.crfDeriv(c, tau1);
		double[] coords = this.getCoordinates().clone();

		coords[tau2] = coords[tau2] + 1E-7;
		NDDOAtom perturbed = new MNDOAtom(this.atomProperties, coords, this.np);

		double newval = perturbed.crfDeriv(c, tau1);

		return 1E7 * (newval - orig);
	}

	@Override
	public double crfAlphaDeriv(NDDOAtom b, int num) {
		MNDOAtom c = (MNDOAtom) b;
		double R = GTO.R(this.getCoordinates(), b.getCoordinates());
		double val = this.atomProperties.getQ() * c.atomProperties.getQ() *
				nom.getG(this.s(), this.s(), b.s(), b.s());

		double returnval = 0;

		if (num == 0 || num == 2) returnval += val * -R / bohr * Math.exp(-this.np.getAlpha() * R / bohr);
		if (num == 1 || num == 2) returnval += val * -R / bohr * Math.exp(-c.np.getAlpha() * R / bohr);

		return returnval;
	}
}
