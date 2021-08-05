package nddoparam;

import nddoparam.mndo.MNDOParams;
import org.jblas.DoubleMatrix;
import runcycle.input.RawAtom;
import runcycle.input.RawMolecule;

import java.util.Arrays;

public abstract class Solution {
	// TODO make most of these private
	public static int maxParamNum = 8;
	public double energy, homo, lumo, hf, dipole;
	public double[] chargedip, hybridip, dipoletot;
	public int charge, mult, nElectrons, nOrbitals;
	public int[][] missingIndices, orbitalIndices;
	public NDDOAtom[] atoms;
	public NDDO6G[] orbitals;
	public int[] atomicNumbers, orbitalAtomNumbers;
	public DoubleMatrix H;
	protected RawMolecule rm;

	public static Solution of(RawMolecule rm, RawAtom[] ras, NDDOParams[] params) {
		NDDOAtom[] atoms;
		if (params instanceof MNDOParams[])
			atoms = RawMolecule.toMNDOAtoms(ras, (MNDOParams[]) params);
		else throw new IllegalStateException("Unidentified params type!");

		if (rm.restricted) return new SolutionR(atoms, rm.atomicNumbers,
				rm.charge, rm.nElectrons, rm.nOrbitals);
		else return new SolutionU(atoms, rm.atomicNumbers,
				rm.charge, rm.mult, rm.nElectrons, rm.nOrbitals);
	}

	protected Solution(NDDOAtom[] atoms, int[] atomicNumbers, int charge,
					   int mult, int nElectrons, int nOrbitals) {
		/*
		 solution give 2 things
		 1. query arrays - atomNumber returns which atom an orbital
		 corresponds to. orbitalIndices returns all orbitals which are from
		 a particular atom. missingIndex returns all orbitals which are not
		 from that atom.

		 2. Fill up H.
		*/
		this.atoms = atoms;
		this.atomicNumbers = atomicNumbers;
		this.charge = charge;
		this.mult = mult;
		this.nElectrons = nElectrons;
		this.nOrbitals = nOrbitals;

		orbitals = new NDDO6G[nOrbitals];
		// give every orbital a unique index
		orbitalIndices = new int[atoms.length][4];
		// orbital indexes of orbitals that an atom doesn't have
		missingIndices = new int[atoms.length][4 * (atoms.length - 1)];
		for (int[] index : missingIndices) Arrays.fill(index, -1);
		// corresponding atom numbers of an orbital
		orbitalAtomNumbers = new int[nOrbitals];

		int overallOrbitalIndex = 0;
		for (int atomIndex = 0, atomsLength = atoms.length;
			 atomIndex < atomsLength; atomIndex++) {
			NDDOAtom atom = atoms[atomIndex];
			int orbitalIndex = 0;
			for (NDDO6G orbital : atom.getOrbitals()) {
				orbitals[overallOrbitalIndex] = orbital;
				orbitalAtomNumbers[overallOrbitalIndex] = atomIndex;
				orbitalIndices[atomIndex][orbitalIndex] = overallOrbitalIndex;
				overallOrbitalIndex++;
				orbitalIndex++;
			}

			if (atom.getAtomProperties().getZ() == 1) {
				orbitalIndices[atomIndex][1] = -1;
				orbitalIndices[atomIndex][2] = -1;
				orbitalIndices[atomIndex][3] = -1;
			}
		}

		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[]{orbitalIndices[j][0],
					orbitalIndices[j][1],
					orbitalIndices[j][2],
					orbitalIndices[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (k != nums[0] && k != nums[1] &&
						k != nums[2] && k != nums[3]) {
					missingIndices[j][counter] = k;
					counter++;
				}
			}
		}

		H = new DoubleMatrix(orbitals.length, orbitals.length);

		// filling up the core matrix in accordance with NDDO formalism
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (k == j) {
					double Huu = orbitals[j].U();

					for (int an = 0; an < atoms.length; an++) {
						if (orbitalAtomNumbers[j] != an) {
							Huu += atoms[an].V(orbitals[j], orbitals[k]);
						}
					}

					H.put(j, k, Huu);
				}
				else if (orbitalAtomNumbers[j] ==
						orbitalAtomNumbers[k]) {
					double Huv = 0;

					for (int an = 0; an < atoms.length; an++) {
						if (orbitalAtomNumbers[j] != an) {
							Huv += atoms[an].V(orbitals[j], orbitals[k]);
						}
					}

					H.put(j, k, Huv);
					H.put(k, j, Huv);
				}
				else {
					double Huk = NDDO6G.beta(orbitals[j], orbitals[k]);

					H.put(j, k, Huk);
					H.put(k, j, Huk);
				}
			}
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

	public abstract DoubleMatrix alphaDensity();

	public abstract DoubleMatrix betaDensity();

	public abstract DoubleMatrix densityMatrix();
}
