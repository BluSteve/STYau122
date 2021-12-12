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
	protected SimpleMatrix H, densityMatrix, alphaDensity, betaDensity;

	protected Solution(MoleculeInfo rm, NDDOAtom[] atoms) {
		this.atoms = atoms;
		this.rm = rm;
		this.charge = rm.charge;
		this.mult = rm.mult;
		this.atomicNumbers = rm.atomicNumbers;
		this.nElectrons = rm.nElectrons;
		this.nOrbitals = rm.nOrbitals;

		orbitals = new NDDO6G[nOrbitals];
		orbsOfAtom = rm.orbsOfAtom; // variable length 2d array
		missingOfAtom = rm.missingOfAtom; // variable length 2d array, complement of orbsOfAtom
		atomOfOrb = rm.atomOfOrb; // nOrbitals length 1d array

		int overallOrbitalIndex = 0;
		for (NDDOAtom atom : atoms) {
			for (NDDO6G orbital : atom.getOrbitals()) {
				orbitals[overallOrbitalIndex] = orbital;
				overallOrbitalIndex++;
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

	protected void findDipole() {
		double[] populations = new double[atoms.length];

		for (int j = 0; j < atoms.length; j++) {
			double sum = 0;
			for (int k : orbsOfAtom[j]) {
				if (k > -1) {
					sum += densityMatrix().get(k, k);
				}
			}

			populations[j] = atoms[j].getAtomProperties().getQ() - sum;
		}

		double[] com = new double[]{0, 0, 0};
		double mass = 0;
		for (NDDOAtom atom : atoms) {
			com[0] += atom.getMass() * atom.getCoordinates()[0];
			com[1] += atom.getMass() * atom.getCoordinates()[1];
			com[2] += atom.getMass() * atom.getCoordinates()[2];
			mass += atom.getMass();
		}

		com[0] /= mass;
		com[1] /= mass;
		com[2] /= mass;

		chargedip = new double[]{0, 0, 0};

		double v = 2.5416;
		for (int j = 0; j < atoms.length; j++) {
			double v1 = v * populations[j];
			chargedip[0] += v1 * (atoms[j].getCoordinates()[0] - com[0]);
			chargedip[1] += v1 * (atoms[j].getCoordinates()[1] - com[1]);
			chargedip[2] += v1 * (atoms[j].getCoordinates()[2] - com[2]);
		}

		hybridip = new double[]{0, 0, 0};

		for (int j = 0; j < atoms.length; j++) {
			if (orbsOfAtom[j].length > 1) { // exclude hydrogen
				double v1 = v * 2 * atoms[j].D1;
				hybridip[0] -= v1 * densityMatrix().get(orbsOfAtom[j][0], orbsOfAtom[j][1]);
				hybridip[1] -= v1 * densityMatrix().get(orbsOfAtom[j][0], orbsOfAtom[j][2]);
				hybridip[2] -= v1 * densityMatrix().get(orbsOfAtom[j][0], orbsOfAtom[j][3]);
			}
		}

		dipoletot = new double[]{chargedip[0] + hybridip[0], chargedip[1] + hybridip[1], chargedip[2] + hybridip[2]};

		dipole = Math.sqrt(dipoletot[0] * dipoletot[0] + dipoletot[1] * dipoletot[1] + dipoletot[2] * dipoletot[2]);
	}

	public MoleculeInfo getRm() {
		return rm;
	}

	public SimpleMatrix densityMatrix() {
		return densityMatrix;
	}

	public SimpleMatrix alphaDensity() {
		return alphaDensity;
	}

	public SimpleMatrix betaDensity() {
		return betaDensity;
	}

	public abstract Solution withNewAtoms(NDDOAtom[] newAtoms);

	public abstract Solution compute();

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
