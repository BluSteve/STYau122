package nddo.solution;

import nddo.NDDO6G;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.mndo.MNDOAtom;
import org.ejml.simple.SimpleMatrix;
import runcycle.input.RawAtom;
import runcycle.input.RawMolecule;
import scf.AtomHandler;

import java.util.Arrays;

public abstract class Solution {
	public static int maxParamNum = 8; // todo compute this on the fly
	private final RawMolecule rm;
	public double energy, homo, lumo, hf, dipole;
	public double[] chargedip, hybridip, dipoletot;
	public int charge, mult, nElectrons, nOrbitals;
	public int[][] missingOfAtom, orbsOfAtom;
	public int[] atomicNumbers, atomOfOrb;
	public NDDOAtom[] atoms;
	public NDDO6G[] orbitals;
	protected SimpleMatrix H;

	protected Solution(RawMolecule rm, NDDOAtom[] atoms) {
		this.atoms = atoms;
		this.rm = rm;
		this.charge = rm.charge;
		this.mult = rm.mult;
		this.atomicNumbers = rm.atomicNumbers;
		this.nElectrons = rm.nElectrons;
		this.nOrbitals = rm.nOrbitals;

		orbitals = new NDDO6G[nOrbitals];
		orbsOfAtom = new int[atoms.length][4];
		atomOfOrb = new int[nOrbitals];
		missingOfAtom = new int[atoms.length][4 * (atoms.length - 1)];
		for (int[] index : missingOfAtom) Arrays.fill(index, -1);

		int overallOrbitalIndex = 0;
		for (int atomIndex = 0, atomsLength = atoms.length;
			 atomIndex < atomsLength; atomIndex++) {
			NDDOAtom atom = atoms[atomIndex];
			int orbitalIndex = 0;
			for (NDDO6G orbital : atom.getOrbitals()) {
				orbitals[overallOrbitalIndex] = orbital;
				atomOfOrb[overallOrbitalIndex] = atomIndex;
				orbsOfAtom[atomIndex][orbitalIndex] = overallOrbitalIndex;
				overallOrbitalIndex++;
				orbitalIndex++;
			}

			if (atom.getAtomProperties().getZ() == 1) {
				orbsOfAtom[atomIndex][1] = -1;
				orbsOfAtom[atomIndex][2] = -1;
				orbsOfAtom[atomIndex][3] = -1;
			}
		}

		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[]{orbsOfAtom[j][0],
					orbsOfAtom[j][1],
					orbsOfAtom[j][2],
					orbsOfAtom[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (k != nums[0] && k != nums[1] &&
						k != nums[2] && k != nums[3]) {
					missingOfAtom[j][counter] = k;
					counter++;
				}
			}
		}

		// filling up the core matrix in accordance with NDDO formalism
		H = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (k == j) {
					double Huu = orbitals[j].U();

					for (int an = 0; an < atoms.length; an++) {
						if (atomOfOrb[j] != an) {
							Huu += atoms[an].V(orbitals[j], orbitals[k]);
						}
					}

					H.set(j, k, Huu);
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					double Huv = 0;

					for (int an = 0; an < atoms.length; an++) {
						if (atomOfOrb[j] != an) {
							Huv += atoms[an].V(orbitals[j], orbitals[k]);
						}
					}

					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else {
					double Huk = NDDO6G.beta(orbitals[j], orbitals[k]);

					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}
	}

	/**
	 * Initializes SCF solution of a molecule.
	 *
	 * @param rm    Contains intrinsic and compulsory information about a molecule. The atoms/ expgeom arrays will
	 *              not be modified.
	 * @param atoms List of NDDO atoms.
	 */
	public static Solution of(RawMolecule rm, NDDOAtom[] atoms) {
		if (rm.restricted) return new SolutionR(rm, atoms).compute();
		else return new SolutionU(rm, atoms).compute();
	}

	public static int[] getNIntegrals(RawMolecule rm) {
		NDDOParams placeholder = new NDDOParams(new double[13]);

		NDDOAtom[] atoms = new MNDOAtom[rm.atoms.length];
		for (int i = 0; i < rm.atoms.length; i++) {
			RawAtom ra = rm.atoms[i];
			atoms[i] = new MNDOAtom(AtomHandler.atoms[ra.Z], ra.coords, placeholder);
		}

		if (rm.restricted) {
			return new int[]{new SolutionR(rm, atoms).findNIntegrals()};
		}
		else {
			SolutionU s = new SolutionU(rm, atoms);
			return new int[]{s.findNCoulombInts(), s.findNExchangeInts()};
		}
	}

	/**
	 * Checks if two DoubleMatrices are similar below a threshold.
	 */
	protected static boolean isSimilar(SimpleMatrix x, SimpleMatrix y,
									   double limit) {
		for (int i = 0; i < y.numRows(); i++) {
			for (int j = 0; j < y.numCols(); j++) {
				if (Math.abs(x.get(i, j) - y.get(i, j)) > limit) {
					return false;
				}
			}
		}

		return true;
	}

	public Solution withNewAtoms(NDDOAtom[] newAtoms) {
		if (this instanceof SolutionR)
			return new SolutionR(rm, newAtoms).compute();
		else if (this instanceof SolutionU)
			return new SolutionU(rm, newAtoms).compute();
		else throw new IllegalStateException("Unidentified Solution type!");
	}

	public RawMolecule getRm() {
		return rm;
	}

	public abstract Solution compute();

	public abstract SimpleMatrix alphaDensity();

	public abstract SimpleMatrix betaDensity();

	public abstract SimpleMatrix densityMatrix();

	@Override
	public String toString() {
		return "Solution{" +
				"rm=" + rm +
				", energy=" + energy +
				", homo=" + homo +
				", lumo=" + lumo +
				", hf=" + hf +
				", dipole=" + dipole +
				", charge=" + charge +
				", mult=" + mult +
				", nElectrons=" + nElectrons +
				", nOrbitals=" + nOrbitals +
				", missingOfAtom=" + Arrays.toString(missingOfAtom) +
				", orbsOfAtom=" + Arrays.toString(orbsOfAtom) +
				", atomicNumbers=" + Arrays.toString(atomicNumbers) +
				", atomOfOrb=" + Arrays.toString(atomOfOrb) +
				", atoms=" + Arrays.toString(atoms) +
				", orbitals=" + Arrays.toString(orbitals) +
				", H=" + H +
				", chargedip=" + Arrays.toString(chargedip) +
				", hybridip=" + Arrays.toString(hybridip) +
				", dipoletot=" + Arrays.toString(dipoletot) +
				'}';
	}
}
