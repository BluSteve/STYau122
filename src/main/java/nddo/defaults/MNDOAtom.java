package nddo.defaults;

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

	private static double getf(MNDOAtom a, MNDOAtom b, double R) {
		return 1 + R / bohr * Pow.exp(-a.np.getAlpha() * R / bohr) +
				Pow.exp(-b.np.getAlpha() * R / bohr);
	}

	private static double getfPrime(MNDOAtom a, MNDOAtom b, double R, int tau) {
		return (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr) *
				Pow.exp(-a.np.getAlpha() * R / bohr) -
				a.np.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (bohr * bohr) *
						Pow.exp(-a.np.getAlpha() * R / bohr) -
				b.np.getAlpha() / bohr * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
						Pow.exp(-b.np.getAlpha() * R / bohr);
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
			f = getf(this, b, R);
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && atomProperties.getZ() == 1)
			f = getf(b, this, R);
		else f = 1 + Pow.exp(-b.np.getAlpha() * R / bohr) + Pow.exp(-this.np.getAlpha() * R / bohr);

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
			f = getf(this, b, R);
			fprime = getfPrime(this, b, R, tau);
		}
		else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && atomProperties.getZ() == 1) {
			f = getf(b, this, R);
			fprime = -getfPrime(b, this, R, tau);
		}
		else {
			f = 1 + Pow.exp(-b.np.getAlpha() * r) + Pow.exp(-this.np.getAlpha() * r);
			fprime = -b.np.getAlpha() / bohr * (coordinates[tau] - b.getCoordinates()[tau]) / R *
					Pow.exp(-b.np.getAlpha() * r) -
					this.np.getAlpha() / bohr * (coordinates[tau] - b.getCoordinates()[tau]) / R *
							Pow.exp(-this.np.getAlpha() * r);
		}

		return fprime * atomProperties.getQ() * b.atomProperties.getQ() * nom.G(this.s(), this.s(), b.s(), b.s()) +
				f * atomProperties.getQ() * b.atomProperties.getQ() * nom.Ggd(this.s(), this.s(), b.s(), b.s(), tau);
	}

	@Override
	public double crfg2d(MNDOAtom c, int tau1, int tau2) {
		double orig = this.crfgd(c, tau1);
		double[] coords = coordinates.clone();

		coords[tau2] = coords[tau2] + 1E-7;
		MNDOAtom perturbed = new MNDOAtom(atomProperties, coords, this.np);

		double newval = perturbed.crfgd(c, tau1);

		return 1E7 * (newval - orig);
	}

	@Override
	public double crfalphapd(MNDOAtom c, int num) {
		double R = GTO.R(coordinates, c.getCoordinates());
		double val = atomProperties.getQ() * c.atomProperties.getQ() * nom.G(this.s(), this.s(), c.s(), c.s());

		double returnval = 0;

		if (num == 0 || num == 2) returnval += val * -R / bohr * Pow.exp(-this.np.getAlpha() * R / bohr);
		if (num == 1 || num == 2) returnval += val * -R / bohr * Pow.exp(-c.np.getAlpha() * R / bohr);

		return returnval;
	}
}
