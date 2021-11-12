package nddoparam;

import org.ejml.simple.SimpleMatrix;
import runcycle.input.RawMolecule;

import java.util.Arrays;

public abstract class Solution {
	// TODO make most of these private
	public static int maxParamNum = 8;
	public double energy, homo, lumo, hf, dipole;
	public double[] chargedip, hybridip, dipoletot;
	public int charge, multiplicity;
	public int[][] missingIndex, orbitalIndices;
	public NDDOAtom[] atoms;
	public int[] atomNumber;
	public double damp = 0.8;
	public int nElectrons;
	public SimpleMatrix H;
	public NDDO6G[] orbitals;
	public int[] atomicNumbers;
	protected RawMolecule rm;

	public Solution(NDDOAtom[] atoms, int charge) {
		/*
		 solution give 2 things
		 1. query arrays - atomNumber returns which atom an orbital
		 corresponds to. orbitalIndices returns all orbitals which are from
		 a particular atom. missingIndex returns all orbitals which are not
		 from that atom.

		 2. Fill up H.
		*/
		this.atoms = atoms;
		this.charge = charge;
		int numOrbitals = 0;

		for (NDDOAtom a : atoms) {
			nElectrons += a.getAtomProperties().getQ();
			numOrbitals += a.getOrbitals().length;
		}
		nElectrons -= charge;

		atomicNumbers = new int[atoms.length];
		for (int num = 0; num < atoms.length; num++) {
			atomicNumbers[num] = atoms[num].getAtomProperties().getZ();
		}

		orbitals = new NDDO6G[numOrbitals];
		orbitalIndices = new int[atoms.length][4];
		atomNumber = new int[orbitals.length];
		int atomIndex = 0;
		int overallIndex = 0;
		for (NDDOAtom atom : atoms) {
			int orbitalIndex = 0;
			for (NDDO6G orbital : atom.getOrbitals()) {
				orbitals[overallIndex] = orbital;
				orbitalIndices[atomIndex][orbitalIndex] = overallIndex;
				atomNumber[overallIndex] = atomIndex;
				overallIndex++;
				orbitalIndex++;
			}

			if (atom.getAtomProperties().getZ() == 1) {
				orbitalIndices[atomIndex][1] = -1;
				orbitalIndices[atomIndex][2] = -1;
				orbitalIndices[atomIndex][3] = -1;
			}
			atomIndex++;
		}

		missingIndex = new int[atoms.length][4 * atoms.length - 4];

		for (int j = 0; j < atoms.length; j++) {
			for (int k = 0; k < 4 * atoms.length - 4; k++) {
				missingIndex[j][k] = -1;
			}
		}

		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[]{orbitalIndices[j][0], orbitalIndices[j][1],
					orbitalIndices[j][2],
					orbitalIndices[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (nums[0] != k && nums[1] != k && nums[2] != k &&
						nums[3] != k) {
					missingIndex[j][counter] = k;
					counter++;
				}
			}
		}

		H = new SimpleMatrix(orbitals.length, orbitals.length);

		//filling up the core matrix in accordance with NDDO formalism

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (k == j) {
					double Huu = orbitals[j].U();

					for (int an = 0; an < atoms.length; an++) {
						if (atomNumber[j] != an) {
							Huu += atoms[an].V(orbitals[j], orbitals[k]);
						}
					}

					H.set(j, k, Huu);
				}
				else if (atomNumber[j] == atomNumber[k]) { // case 2
					double Huv = 0;
					for (int an = 0; an < atoms.length; an++) {
						if (atomNumber[j] != an) {
							Huv += atoms[an].V(orbitals[j], orbitals[k]);
						}
					}
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { // case 3
					double Huk = NDDO6G.beta(orbitals[j], orbitals[k]);
					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}
	}

	/**
	 * Checks if two DoubleMatrices are similar below a threshold.
	 * @param x
	 * @param y
	 * @param limit
	 * @return
	 */
	public static boolean isSimilar(SimpleMatrix x, SimpleMatrix y,
									double limit) {
		for (int i = 0; i < y.numRows(); i++) {
			for (int j = 0; j < y.numCols(); j++) {
				if (Math.abs(x.get(i, j) - y.get(i, j)) > limit) {
					System.err.println (i + ", " + j + ": " + Math.abs(x.get(i, j) - y.get(i, j)));
					return false;
				}
			}
		}
		return true;
	}

	public RawMolecule getRm() {
		return rm;
	}

	public abstract Solution setRm(RawMolecule rm);

	public abstract SimpleMatrix alphaDensity();

	public abstract SimpleMatrix betaDensity();

	public abstract SimpleMatrix densityMatrix();

	@Override
	public String toString() {
		return "Solution{" +
				", energy=" + energy +
				", homo=" + homo +
				", lumo=" + lumo +
				", hf=" + hf +
				", dipole=" + dipole +
				", chargedip=" + Arrays.toString(chargedip) +
				", hybridip=" + Arrays.toString(hybridip) +
				", dipoletot=" + Arrays.toString(dipoletot) +
				", charge=" + charge +
				", multiplicity=" + multiplicity +
				", missingIndex=" + Arrays.toString(missingIndex) +
				", index=" + Arrays.toString(orbitalIndices) +
				", atoms=" + Arrays.toString(atoms) +
				", atomNumber=" + Arrays.toString(atomNumber) +
				", damp=" + damp +
				", nElectrons=" + nElectrons +
				", H=" + H +
				", orbitals=" + Arrays.toString(orbitals) +
				", atomicNumbers=" + Arrays.toString(atomicNumbers) +
				'}';
	}
}
