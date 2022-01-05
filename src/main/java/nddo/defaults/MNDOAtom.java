package nddo.defaults;

import nddo.Constants;
import nddo.NDDOParams;
import nddo.scf.GTO;
import nddo.structs.AtomProperties;
import nddo.structs.OrbitalProperties;
import tools.Pow;

import static nddo.Constants.bohr;
import static nddo.State.nom;

public class MNDOAtom extends NDDOAtomBasic<MNDOAtom> {
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

	private static double getfhydrogen(MNDOAtom a, MNDOAtom b, double R) {
		return 1 + R / bohr * Pow.exp(-a.np.getAlpha() * R / bohr) +
				Pow.exp(-b.np.getAlpha() * R / bohr);
	}

	private static double getf(MNDOAtom a, MNDOAtom b, double R) {
		return 1 + Pow.exp(-b.np.getAlpha() * R / bohr) + Pow.exp(-a.np.getAlpha() * R / bohr);
	}

	private static double getfhydrogenpd(MNDOAtom a, double R) {
		return -R * R / (bohr * bohr) * Pow.exp(-a.np.getAlpha() * R / bohr);
	}

	private static double getfhydrogenp2d(MNDOAtom a, double R) {
		return R * R * R / (bohr * bohr * bohr) * Pow.exp(-a.np.getAlpha() * R / bohr);
	}

	private static double getfpd(MNDOAtom a, double R) {
		return -R / bohr * Pow.exp(-a.np.getAlpha() * R / bohr);
	}

	private static double getfp2d(MNDOAtom a, double R) {
		return R * R / (bohr * bohr) * Pow.exp(-a.np.getAlpha() * R / bohr);
	}

	private static double getfhydrogenpgd(MNDOAtom a, MNDOAtom b, double R, int tau) {
		return -2 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (bohr * bohr) * Pow.exp(-a.np.getAlpha() * R / bohr)
				+ a.np.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) * R / (bohr * bohr * bohr) *
				Pow.exp(-a.np.getAlpha() * R / bohr);
	}

	private static double getfhydrogenp2gd(MNDOAtom a, MNDOAtom b, double R, int tau) {
		return 3 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) * R / (bohr * bohr * bohr) *
				Pow.exp(-a.np.getAlpha() * R / bohr)
				- a.np.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) * R * R / (bohr * bohr * bohr * bohr) *
				Pow.exp(-a.np.getAlpha() * R / bohr);
	}

	private static double getfpgd(MNDOAtom a, MNDOAtom c, double R, int tau) {
		return -(a.getCoordinates()[tau] - c.getCoordinates()[tau]) / (R * bohr) * Pow.exp(-a.np.getAlpha() * R / bohr)
				+ a.np.getAlpha() * (a.getCoordinates()[tau] - c.getCoordinates()[tau]) / (bohr * bohr) *
				Pow.exp(-a.np.getAlpha() * R / bohr);
	}

	private static double getfhydrogengd(MNDOAtom a, MNDOAtom b, double R, int tau) {
		return (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr) * Pow.exp(-a.np.getAlpha() * R / bohr)
				- a.np.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (bohr * bohr) *
				Pow.exp(-a.np.getAlpha() * R / bohr)
				- b.np.getAlpha() / bohr * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
				Pow.exp(-b.np.getAlpha() * R / bohr);
	}

	private static double getfgd(MNDOAtom a, MNDOAtom b, double R, int tau) {
		return -b.np.getAlpha() / bohr * (a.coordinates[tau] - b.getCoordinates()[tau]) / R *
				Pow.exp(-b.np.getAlpha() * R / bohr) -
				a.np.getAlpha() / bohr * (a.coordinates[tau] - b.getCoordinates()[tau]) / R *
						Pow.exp(-a.np.getAlpha() * R / bohr);
	}

	private static double getfp2gd(MNDOAtom a, MNDOAtom c, double R, int tau) {
		return 2 * (a.getCoordinates()[tau] - c.getCoordinates()[tau]) / (bohr * bohr) *
				Pow.exp(-a.np.getAlpha() * R / bohr)
				- a.np.getAlpha() * (a.getCoordinates()[tau] - c.getCoordinates()[tau]) * R / (bohr * bohr * bohr) *
				Pow.exp(-a.np.getAlpha() * R / bohr);
	}



	@Override
	public MNDOAtom withNewParams(NDDOParams np) {
		return new MNDOAtom(atomProperties, coordinates, np);
	}

	@Override
	public MNDOAtom withNewCoords(double[] coordinates) {
		return new MNDOAtom(atomProperties, coordinates, np);
	}

	@Override
	public MNDOAtom copy() {
		return new MNDOAtom(atomProperties, coordinates, np);
	}

	@Override
	public double crf(MNDOAtom b) {
		double f;
		double R = GTO.R(coordinates, b.getCoordinates());

		if ((atomProperties.getZ() == 7 || atomProperties.getZ() == 8) && b.atomProperties.getZ() == 1)
			f = getfhydrogen(this, b, R);
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && atomProperties.getZ() == 1)
			f = getfhydrogen(b, this, R);
		else f = getf(this, b, R);

		return f * atomProperties.getQ() * b.atomProperties.getQ() * nom.G(this.s(), this.s(), b.s(),
				b.s());
	}

	@Override
	public double crfgd(MNDOAtom b, int tau) {
		double f;
		double fprime;
		double R = GTO.R(coordinates, b.getCoordinates());
		double r = R / bohr;

		if ((atomProperties.getZ() == 7 || atomProperties.getZ() == 8) && b.atomProperties.getZ() == 1) {
			f = getfhydrogen(this, b, R);
			fprime = getfhydrogengd(this, b, R, tau);
		}
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && atomProperties.getZ() == 1) {
			f = getfhydrogen(b, this, R);
			fprime = -getfhydrogengd(b, this, R, tau);
		}
		else {
			f = 1 + Pow.exp(-b.np.getAlpha() * r) + Pow.exp(-this.np.getAlpha() * r);
			fprime = getfgd(this, b, R, tau);
		}

		return fprime * atomProperties.getQ() * b.atomProperties.getQ() * nom.G(this.s(), this.s(), b.s(), b.s()) +
				f * atomProperties.getQ() * b.atomProperties.getQ() * nom.Ggd(this.s(), this.s(), b.s(), b.s(), tau);
	}

	@Override
	public double crfg2d(MNDOAtom c, int tau1, int tau2) {
		double orig = this.crfgd(c, tau1);
		double[] coords = coordinates.clone();

		coords[tau2] += Constants.LAMBDA;
		MNDOAtom perturbed = new MNDOAtom(atomProperties, coords, this.np);

		double newval = perturbed.crfgd(c, tau1);

		return 1 / Constants.LAMBDA * (newval - orig);
	}

	public double crfalphapd(MNDOAtom c, int num) {
		double R = GTO.R(coordinates, c.getCoordinates());
		double val = atomProperties.getQ() * c.atomProperties.getQ() * nom.G(this.s(), this.s(), c.s(), c.s());

		double returnval = 0;


		if (atomProperties.getZ() == 1 && (c.getAtomProperties().getZ() == 7 || c.getAtomProperties().getZ() == 8)) {
			if (num == 0) returnval += val * getfpd(this, R);
			else if (num == 1) returnval += val * getfhydrogenpd(c, R);
		}

		else if ((atomProperties.getZ() == 7 || atomProperties.getZ() == 8) && c.getAtomProperties().getZ() == 1) {
			if (num == 0) returnval += val * getfhydrogenpd(this, R);
			else if (num == 1) returnval += val * getfpd(c, R);
		}

		else {
			if (num == 0 || num == 2) returnval += val * getfpd(this, R);
			if (num == 1 || num == 2) returnval += val * getfpd(c, R);
		}

		return returnval;
	}

	public double crfalphap2d(MNDOAtom c, int num1) {

		double R = GTO.R(coordinates, c.getCoordinates());
		double val = atomProperties.getQ() * c.atomProperties.getQ() * nom.G(this.s(), this.s(), c.s(), c.s());

		double returnval = 0;


		if (atomProperties.getZ() == 1 && (c.getAtomProperties().getZ() == 7 || c.getAtomProperties().getZ() == 8)) {
			if (num1 == 0) returnval += val * getfp2d(this, R);
			else if (num1 == 1) returnval += val * getfhydrogenp2d(c, R);
		}

		else if ((atomProperties.getZ() == 7 || atomProperties.getZ() == 8) && c.getAtomProperties().getZ() == 1) {
			if (num1 == 0) returnval += val * getfhydrogenp2d(this, R);
			else if (num1 == 1) returnval += val * getfp2d(c, R);
		}

		else {
			if (num1 == 0 || num1 == 2) returnval += val * getfp2d(this, R);
			if (num1 == 1 || num1 == 2) returnval += val * getfp2d(c, R);
		}

		return returnval;
	}

	public double crfalphapgd(MNDOAtom c, int num, int tau) {
		double R = GTO.R(coordinates, c.getCoordinates());
		double val = atomProperties.getQ() * c.atomProperties.getQ() * nom.G(this.s(), this.s(), c.s(), c.s());
		double valgd = atomProperties.getQ() * c.atomProperties.getQ() * nom.Ggd(this.s(), this.s(), c.s(), c.s(),
				tau);

		double returnval = 0;

		if (atomProperties.getZ() == 1 && (c.getAtomProperties().getZ() == 7 || c.getAtomProperties().getZ() == 8)) {
			if (num == 0) {
				returnval += valgd * getfpd(this, R);
				returnval += val * getfpgd(this, c, R, tau);
			}
			if (num == 1) {
				returnval += valgd * getfhydrogenpd(c, R);
				returnval += val * -getfhydrogenpgd(c, this, R, tau);
			}
		}

		else if ((atomProperties.getZ() == 7 || atomProperties.getZ() == 8) && c.getAtomProperties().getZ() == 1) {
			if (num == 0) {
				returnval += valgd * getfhydrogenpd(this, R);
				returnval += val * getfhydrogenpgd(this, c, R, tau);
			}
			if (num == 1) {
				returnval += valgd * getfpd(c, R);
				returnval += val * -getfpgd(c, this, R, tau);
			}
		}

		else {
			if (num == 0 || num == 2) {
				returnval += valgd * getfpd(this, R);
				returnval += val * getfpgd(this, c, R, tau);
			}
			if (num == 1 || num == 2) {
				returnval += valgd * getfpd(c, R);
				returnval += val * -getfpgd(c, this, R, tau);
			}
		}

		return returnval;
	}

	@Override
	public double crfalphap2gd(MNDOAtom c, int num, int tau) {
		double R = GTO.R(coordinates, c.getCoordinates());
		double val = atomProperties.getQ() * c.atomProperties.getQ() * nom.G(this.s(), this.s(), c.s(), c.s());
		double valgd = atomProperties.getQ() * c.atomProperties.getQ() * nom.Ggd(this.s(), this.s(), c.s(), c.s(),
				tau);

		double returnval = 0;

		if (atomProperties.getZ() == 1 && (c.getAtomProperties().getZ() == 7 || c.getAtomProperties().getZ() == 8)) {
			if (num == 0) {
				returnval += valgd * getfp2d(this, R);
				returnval += val * getfp2gd(this, c, R, tau);
			}
			if (num == 1) {
				returnval += valgd * getfhydrogenp2d(c, R);
				returnval += val * -getfhydrogenp2gd(c, this, R, tau);
			}
		}

		else if ((atomProperties.getZ() == 7 || atomProperties.getZ() == 8) && c.getAtomProperties().getZ() == 1) {
			if (num == 0) {
				returnval += valgd * getfhydrogenp2d(this, R);
				returnval += val * getfhydrogenp2gd(this, c, R, tau);
			}
			if (num == 1) {
				returnval += valgd * getfp2d(c, R);
				returnval += val * -getfp2gd(c, this, R, tau);
			}
		}

		else {
			if (num == 0 || num == 2) {
				returnval += valgd * getfpd(this, R);
				returnval += val * getfpgd(this, c, R, tau);
			}
			if (num == 1 || num == 2) {
				returnval += valgd * getfpd(c, R);
				returnval += val * -getfpgd(c, this, R, tau);
			}
		}

		return returnval;
	}
}
