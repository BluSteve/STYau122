package nddo.solution;

import nddo.MoleculeInfo;
import nddo.NDDO6G;
import nddo.NDDOAtom;
import org.ejml.simple.SimpleMatrix;

import java.util.Arrays;

public abstract class Solution {
	public static int maxParamNum = 8; // todo compute this on the fly
	protected final MoleculeInfo rm;
	public double energy, homo, lumo, hf, dipole;
	public double[] chargedip, hybridip, dipoletot;
	public int charge, mult, nElectrons, nOrbitals;
	public int[][] missingOfAtom, orbsOfAtom;
	public int[] atomicNumbers, atomOfOrb;
	public NDDOAtom[] atoms;
	public NDDO6G[] orbitals;
	protected SimpleMatrix H;


	protected Solution(MoleculeInfo rm, NDDOAtom[] atoms) {
		this.atoms = atoms;
		this.rm = rm;
		this.charge = rm.charge;
		this.mult = rm.mult;
		this.atomicNumbers = rm.atomicNumbers;
		this.nElectrons = rm.nElectrons;
		this.nOrbitals = rm.nOrbitals;

		orbitals = new NDDO6G[nOrbitals];
		orbsOfAtom = rm.orbsOfAtom;
		atomOfOrb = rm.atomOfOrb;
		missingOfAtom = rm.missingOfAtom;

		int overallOrbitalIndex = 0;
		for (NDDOAtom atom : atoms) {
			for (NDDO6G orbital : atom.getOrbitals()) {
				orbitals[overallOrbitalIndex] = orbital;
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
	 * @param mi    Contains intrinsic and compulsory information about a molecule.
	 * @param atoms List of NDDO atoms.
	 */
	public static Solution of(MoleculeInfo mi, NDDOAtom[] atoms) {
		if (mi.restricted) return new SolutionR(mi, atoms).compute();
		else return new SolutionU(mi, atoms).compute();
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

	public MoleculeInfo getRm() {
		return rm;
	}

	public abstract Solution withNewAtoms(NDDOAtom[] newAtoms);

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
