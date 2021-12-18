package nddo.defaults;

import nddo.NDDOParams;
import nddo.scf.GTO;
import nddo.structs.AtomProperties;
import nddo.structs.OrbitalProperties;

import static nddo.Constants.bohr;
import static nddo.State.nom;

public class MNDOAtom extends NDDOAtomBasic {
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

	private static double getf(NDDOAtomBasic a, NDDOAtomBasic b, double R) {
		return 1 + R / bohr * Math.exp(-a.getParams().getAlpha() * R / bohr) +
				Math.exp(-b.getParams().getAlpha() * R / bohr);
	}

	private static double getfPrime(NDDOAtomBasic a, NDDOAtomBasic b, double R, int tau) {
		return (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * bohr) *
				Math.exp(-a.getParams().getAlpha() * R / bohr) -
				a.getParams().getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (bohr * bohr) *
						Math.exp(-a.getParams().getAlpha() * R / bohr) -
				b.getParams().getAlpha() / bohr * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R *
						Math.exp(-b.getParams().getAlpha() * R / bohr);
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
	public double crf(NDDOAtomBasic b) {
		double f;
		double R = GTO.R(coordinates, b.getCoordinates());

		if ((atomProperties.getZ() == 7 || atomProperties.getZ() == 8) && b.getAtomProperties().getZ() == 1)
			f = getf(this, b, R);
		else if ((b.getAtomProperties().getZ() == 7 || b.getAtomProperties().getZ() == 8) && atomProperties.getZ() == 1)
			f = getf(b, this, R);
		else f = 1 + Math.exp(-b.getParams().getAlpha() * R / bohr) + Math.exp(-this.np.getAlpha() * R / bohr);

		return f * atomProperties.getQ() * b.getAtomProperties().getQ() * nom.G(this.s(), this.s(), b.s(),
				b.s());
	}

	@Override
	public double crfgd(NDDOAtomBasic b, int tau) {
		double f;
		double fprime;
		double R = GTO.R(coordinates, b.getCoordinates());
		double r = R / bohr;

		if ((atomProperties.getZ() == 7 || atomProperties.getZ() == 8) && b.getAtomProperties().getZ() == 1) {
			f = getf(this, b, R);
			fprime = getfPrime(this, b, R, tau);
		}
		else if ((b.getAtomProperties().getZ() == 7 || b.getAtomProperties().getZ() == 8) &&
				atomProperties.getZ() == 1) {
			f = getf(b, this, R);
			fprime = -getfPrime(b, this, R, tau);
		}
		else {
			f = 1 + Math.exp(-b.getParams().getAlpha() * r) + Math.exp(-this.np.getAlpha() * r);
			fprime = -b.getParams().getAlpha() / bohr * (coordinates[tau] - b.getCoordinates()[tau]) / R *
					Math.exp(-b.getParams().getAlpha() * r) -
					this.np.getAlpha() / bohr * (coordinates[tau] - b.getCoordinates()[tau]) / R *
							Math.exp(-this.np.getAlpha() * r);
		}

		return fprime * atomProperties.getQ() * b.getAtomProperties().getQ() *
				nom.G(this.s(), this.s(), b.s(), b.s()) +
				f * atomProperties.getQ() * b.getAtomProperties().getQ() *
						nom.Ggd(this.s(), this.s(), b.s(), b.s(), tau);
	}

	@Override
	public double crfg2d(NDDOAtomBasic c, int tau1, int tau2) {/*TODO*/
		double orig = this.crfgd(c, tau1); // todo maybe cast instead of using nddoatombasic directly
		double[] coords = coordinates.clone();

		coords[tau2] = coords[tau2] + 1E-7;
		MNDOAtom perturbed = new MNDOAtom(atomProperties, coords, this.np);

		double newval = perturbed.crfgd(c, tau1);

		return 1E7 * (newval - orig);
	}

	@Override
	public double crfpd(NDDOAtomBasic c, int num) {
		double R = GTO.R(coordinates, c.getCoordinates());
		double val = atomProperties.getQ() * c.getAtomProperties().getQ() *
				nom.G(this.s(), this.s(), c.s(), c.s());

		double returnval = 0;

		if (num == 0 || num == 2) returnval += val * -R / bohr * Math.exp(-this.np.getAlpha() * R / bohr);
		if (num == 1 || num == 2) returnval += val * -R / bohr * Math.exp(-c.getParams().getAlpha() * R / bohr);

		return returnval;
	}
}
