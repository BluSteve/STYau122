package nddoparam;

import org.jblas.DoubleMatrix;
import runcycle.input.RawMolecule;

import java.util.Arrays;

public abstract class Solution implements Cloneable {
	// TODO make most of these private
	public static int maxParamNum = 8;
	public double energy, homo, lumo, hf, dipole;
	public double[] chargedip, hybridip, dipoletot;
	public int charge, multiplicity;
	public int[][] missingIndex, index;
	public NDDOAtom[] atoms;
	public int[] atomNumber;
	public double damp = 0.8;
	public int nElectrons;
	public DoubleMatrix H;
	public NDDO6G[] orbitals;
	public int[] atomicNumbers;
	protected RawMolecule rm;

	public Solution(NDDOAtom[] atoms, int charge) {
		for (NDDOAtom a : atoms) {
			nElectrons += a.getAtomProperties().getQ();
		}
		this.atoms = atoms;
		atomicNumbers = new int[atoms.length];
		nElectrons -= charge;
		this.charge = charge;
		int i = 0;

		for (NDDOAtom a : atoms) {
			i += a.getOrbitals().length;
		}

		for (int num = 0; num < atoms.length; num++) {
			atomicNumbers[num] = atoms[num].getAtomProperties().getZ();
		}

		orbitals = new NDDO6G[i];

		i = 0;

		index = new int[atoms.length][4];
		atomNumber = new int[orbitals.length];
		int count = 0;
		int count2;
		for (NDDOAtom a : atoms) {
			count2 = 0;
			for (NDDO6G orbital : a.getOrbitals()) {
				orbitals[i] = orbital;
				index[count][count2] = i;
				atomNumber[i] = count;
				i++;
				count2++;
			}


			if (a.getAtomProperties().getZ() == 1) {
				index[count][1] = -1;
				index[count][2] = -1;
				index[count][3] = -1;
			}
			count++;
		}

		missingIndex = new int[atoms.length][4 * atoms.length - 4];

		for (int j = 0; j < atoms.length; j++) {
			for (int k = 0; k < 4 * atoms.length - 4; k++) {
				missingIndex[j][k] = -1;
			}
		}

		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[]{index[j][0], index[j][1], index[j][2],
					index[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (nums[0] != k && nums[1] != k && nums[2] != k &&
						nums[3] != k) {
					missingIndex[j][counter] = k;
					counter++;
				}
			}
		}

		H = new DoubleMatrix(orbitals.length, orbitals.length);

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

					H.put(j, k, Huu);
				}
				else if (atomNumber[j] == atomNumber[k]) { // case 2
					double Huv = 0;
					for (int an = 0; an < atoms.length; an++) {
						if (atomNumber[j] != an) {
							Huv += atoms[an].V(orbitals[j], orbitals[k]);
						}
					}
					H.put(j, k, Huv);
					H.put(k, j, Huv);
				}
				else { // case 3
					double Huk = NDDO6G.beta(orbitals[j], orbitals[k]);
					H.put(j, k, Huk);
					H.put(k, j, Huk);
				}
			}
		}
	}

	public static boolean isSimilar(DoubleMatrix x, DoubleMatrix y,
									double limit) {
		for (int i = 0; i < y.rows; i++) {
			for (int j = 0; j < y.columns; j++) {
				if (Math.abs(x.get(i, j) - y.get(i, j)) > limit) {
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

	public abstract Solution clone();

	public abstract DoubleMatrix alphaDensity();

	public abstract DoubleMatrix betaDensity();

	public abstract DoubleMatrix densityMatrix();


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
				", index=" + Arrays.toString(index) +
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
