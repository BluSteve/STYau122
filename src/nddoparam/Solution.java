package nddoparam;

import nddoparam.mndo.MNDOParams;
import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;
import runcycle.input.RawAtom;
import runcycle.input.RawMolecule;
import scf.Utils;

import java.util.Arrays;

public abstract class Solution {
	// TODO make most of these private
	public static int maxParamNum = 8;
	public double energy, homo, lumo, hf, dipole;
	public double[] chargedip, hybridip, dipoletot;
	public int charge, mult, nElectrons, nOrbitals;
	public int[][] missingOfAtom, orbsOfAtom;
	public NDDOAtom[] atoms;
	public NDDO6G[] orbitals;
	public int[] atomicNumbers, atomOfOrb;
	public SimpleMatrix H;
	protected RawMolecule rm;

	protected Solution(NDDOAtom[] atoms, RawMolecule rm) {
		this.atoms = atoms;
		this.rm = rm;
		this.atomicNumbers = rm.atomicNumbers;
		this.charge = rm.charge;
		this.mult = rm.mult;
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

		fillH();
	}

	public static Solution of(RawMolecule rm, RawAtom[] ras,
							  NDDOParams[] params) {
		NDDOAtom[] atoms;
		if (params instanceof MNDOParams[])
			atoms = RawMolecule.toMNDOAtoms(ras, (MNDOParams[]) params);
		else throw new IllegalStateException("Unidentified params type!");

		if (rm.restricted) return new SolutionR(atoms, rm).compute();
		else return new SolutionU(atoms, rm).compute();
	}

	public static int[] getNIntegrals(RawMolecule rm) {
		MNDOParams[] placeholder = new MNDOParams[Utils.maxAtomNum];
		for (int i = 0; i < Utils.maxAtomNum; i++) {
			placeholder[i] = new MNDOParams();
		}
		NDDOAtom[] atoms = RawMolecule.toMNDOAtoms(rm.atoms, placeholder);
		if (rm.restricted)
			return new int[]{new SolutionR(atoms, rm).findNIntegrals()};
		else {
			SolutionU s = new SolutionU(atoms, rm);
			return new int[]{s.findNCoulombInts(), s.findNExchangeInts()};
		}
	}

	/**
	 * Checks if two DoubleMatrices are similar below a threshold.
	 *
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
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * filling up the core matrix in accordance with NDDO formalism
	 */
	protected void fillH() {
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

	public Solution withNewAtoms(NDDOAtom[] newAtoms) {
		if (this instanceof SolutionR)
			return new SolutionR(newAtoms, rm).compute();
		else if (this instanceof SolutionU)
			return new SolutionU(newAtoms, rm).compute();
		else throw new IllegalStateException("Unidentified Solution type!");
	}

	public abstract Solution compute();

	public RawMolecule getRm() {
		return rm;
	}

	public abstract DoubleMatrix alphaDensity();

	public abstract DoubleMatrix betaDensity();

	public abstract DoubleMatrix densityMatrix();

	@Override
	public String toString() {
		return "Solution{" +
				"energy=" + energy +
				", homo=" + homo +
				", lumo=" + lumo +
				", hf=" + hf +
				", dipole=" + dipole +
				", chargedip=" + Arrays.toString(chargedip) +
				", hybridip=" + Arrays.toString(hybridip) +
				", dipoletot=" + Arrays.toString(dipoletot) +
				", charge=" + charge +
				", mult=" + mult +
				", nElectrons=" + nElectrons +
				", nOrbitals=" + nOrbitals +
				", missingOfAtom=" + Arrays.toString(missingOfAtom) +
				", orbsOfAtom=" + Arrays.toString(orbsOfAtom) +
				", atoms=" + Arrays.toString(atoms) +
				", orbitals=" + Arrays.toString(orbitals) +
				", atomicNumbers=" + Arrays.toString(atomicNumbers) +
				", atomOfOrb=" + Arrays.toString(atomOfOrb) +
				", H=" + H +
				", rm=" + rm +
				'}';
	}
}
