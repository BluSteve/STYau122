package nddo.solution;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.defaults.NDDO6G;
import nddo.structs.MoleculeInfo;
import org.ejml.simple.SimpleMatrix;
import tools.Pow;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;

import static nddo.State.nom;

public abstract class Solution {
	public final static int maxParamNum = 8; // todo compute this on the fly
	protected	static final int[][][] TBRS = {
			{{}},
			{{0}},
			{{1}, {0, 1}},
			{{2}, {0, 2}, {1, 2}, {0, 1, 2}},
			{{3}, {0, 3}, {1, 3}, {2, 3}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}, {0, 1, 2, 3}},
			{{4}, {0, 4}, {1, 4}, {2, 4}, {3, 4}, {0, 1, 4}, {0, 2, 4}, {0, 3, 4}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4},
					{0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}, {1, 2, 3, 4}, {0, 1, 2, 3, 4}},
			{{5}, {0, 5}, {1, 5}, {2, 5}, {3, 5}, {4, 5}, {0, 1, 5}, {0, 2, 5}, {0, 3, 5}, {0, 4, 5}, {1, 2, 5},
					{1, 3, 5}, {1, 4, 5}, {2, 3, 5}, {2, 4, 5}, {3, 4, 5}, {0, 1, 2, 5}, {0, 1, 3, 5}, {0, 1, 4, 5},
					{0, 2, 3, 5}, {0, 2, 4, 5}, {0, 3, 4, 5}, {1, 2, 3, 5}, {1, 2, 4, 5}, {1, 3, 4, 5}, {2, 3, 4, 5},
					{0, 1, 2, 3, 5}, {0, 1, 2, 4, 5}, {0, 1, 3, 4, 5}, {0, 2, 3, 4, 5}, {1, 2, 3, 4, 5},
					{0, 1, 2, 3, 4, 5}},
			{{6}, {0, 6}, {1, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6}, {0, 1, 6}, {0, 2, 6}, {0, 3, 6}, {0, 4, 6},
					{0, 5, 6}, {1, 2, 6}, {1, 3, 6}, {1, 4, 6}, {1, 5, 6}, {2, 3, 6}, {2, 4, 6}, {2, 5, 6}, {3, 4, 6},
					{3, 5, 6}, {4, 5, 6}, {0, 1, 2, 6}, {0, 1, 3, 6}, {0, 1, 4, 6}, {0, 1, 5, 6}, {0, 2, 3, 6},
					{0, 2, 4, 6}, {0, 2, 5, 6}, {0, 3, 4, 6}, {0, 3, 5, 6}, {0, 4, 5, 6}, {1, 2, 3, 6}, {1, 2, 4, 6},
					{1, 2, 5, 6}, {1, 3, 4, 6}, {1, 3, 5, 6}, {1, 4, 5, 6}, {2, 3, 4, 6}, {2, 3, 5, 6}, {2, 4, 5, 6},
					{3, 4, 5, 6}, {0, 1, 2, 3, 6}, {0, 1, 2, 4, 6}, {0, 1, 2, 5, 6}, {0, 1, 3, 4, 6}, {0, 1, 3, 5, 6},
					{0, 1, 4, 5, 6}, {0, 2, 3, 4, 6}, {0, 2, 3, 5, 6}, {0, 2, 4, 5, 6}, {0, 3, 4, 5, 6},
					{1, 2, 3, 4, 6}, {1, 2, 3, 5, 6}, {1, 2, 4, 5, 6}, {1, 3, 4, 5, 6}, {2, 3, 4, 5, 6},
					{0, 1, 2, 3, 4, 6}, {0, 1, 2, 3, 5, 6}, {0, 1, 2, 4, 5, 6}, {0, 1, 3, 4, 5, 6}, {0, 2, 3, 4, 5, 6},
					{1, 2, 3, 4, 5, 6}, {0, 1, 2, 3, 4, 5, 6}}
	};
	protected static final double[] ZEROS = new double[8];

	public final int charge, mult, nElectrons, nOrbitals;
	public final int[][] missingOfAtom, orbsOfAtom;
	public final int[] atomicNumbers, atomOfOrb;
	public final NDDOAtom[] atoms;
	public final NDDOOrbital[] orbitals;
	public final MoleculeInfo rm;
	protected final SimpleMatrix H;
	public double energy, homo, lumo, hf, dipole;
	public double[] chargedip, hybridip, dipoletot;
	protected SimpleMatrix densityMatrix, alphaDensity, betaDensity;

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
			for (NDDOOrbital orbital : atom.getOrbitals()) {
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
					double Huk = nom.H(orbitals[j], orbitals[k]);

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

	protected static SimpleMatrix commutator(SimpleMatrix F, SimpleMatrix D) {
		return F.mult(D).minusi(D.mult(F));
	}

	protected static SimpleMatrix removeElements(SimpleMatrix original, int[] tbr) {
		int newN = original.numRows() - tbr.length;
		boolean isMatrix = original.numCols() > 1;
		SimpleMatrix newarray = isMatrix ? new SimpleMatrix(newN, newN) : new SimpleMatrix(newN, 1);

		int[] tbk = new int[newN];
		int q = 0;
		for (int i = 0; i < original.numRows(); i++) {
			boolean inside = false;
			for (int index : tbr) {
				if (i == index) {
					inside = true;
					break;
				}
			}

			if (!inside) {
				tbk[q] = i;
				q++;
			}
		}

		if (isMatrix) {
			for (int i = 0; i < newarray.numRows(); i++) {
				for (int j = 0; j < newarray.numCols(); j++) {
					newarray.set(i, j, original.get(tbk[i], tbk[j]));
				}
			}
		}
		else {
			for (int i = 0; i < newarray.numRows(); i++) {
				newarray.set(i, 0, original.get(tbk[i], 0));
			}
		}

		return newarray;
	}

	protected static SimpleMatrix addRows(SimpleMatrix original, int[] tbr) {
		// add zero row at tbr
		ArrayList<Double> array = new ArrayList<>();

		for (double i : original.getDDRM().data) {
			array.add(i);
		}

		for (int index : tbr) {
			array.add(index, 0.0);
		}

		return new SimpleMatrix(original.numRows() + tbr.length, 1,
				true, Utils.toDoubles(array));
	}

	protected static double sigmoid(double x) {
		return 1e-10 / (1 + Pow.exp(-0.1 * (x - 53)));
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
			com[0] += atom.getAtomProperties().getMass() * atom.getCoordinates()[0];
			com[1] += atom.getAtomProperties().getMass() * atom.getCoordinates()[1];
			com[2] += atom.getAtomProperties().getMass() * atom.getCoordinates()[2];
			mass += atom.getAtomProperties().getMass();
		}

		com[0] /= mass;
		com[1] /= mass;
		com[2] /= mass;

		chargedip = new double[]{0, 0, 0};

		for (int j = 0; j < atoms.length; j++) {
			double v1 = Constants.DIPOLECONV * populations[j];
			chargedip[0] += v1 * (atoms[j].getCoordinates()[0] - com[0]);
			chargedip[1] += v1 * (atoms[j].getCoordinates()[1] - com[1]);
			chargedip[2] += v1 * (atoms[j].getCoordinates()[2] - com[2]);
		}

		hybridip = new double[]{0, 0, 0};

		for (int j = 0; j < atoms.length; j++) {
			if (orbsOfAtom[j].length > 1) { // exclude hydrogen
				double v1 = Constants.DIPOLECONV * 2 * atoms[j].D1();
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
