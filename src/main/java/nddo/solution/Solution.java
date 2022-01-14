package nddo.solution;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.defaults.NDDO6G;
import nddo.structs.MoleculeInfo;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;
import tools.Pow;
import tools.Utils;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import static nddo.State.config;
import static nddo.State.nom;
import static org.ejml.dense.row.CommonOps_DDRM.elementMin;

public abstract class Solution {
	public final static int maxParamNum = 8; // todo compute this on the fly
	protected static final int[][][] TBRS = {
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
	protected static final int LENGTH = 8;
	protected static final int LENGTH1 = LENGTH - 1;
	protected static final double[] ZEROS = new double[LENGTH];
	protected static final double[] THRESHOLDS =
			new double[10 * Math.max(config.rhf_numIt_max, config.uhf_numIt_max) + 1];

	static {
		for (int i = 0; i < THRESHOLDS.length; i++) {
			THRESHOLDS[i] = 1e-10 / (1 + Pow.exp(-0.1 * (i - 53)));
		}
	}

	public final int charge, mult, nElectrons, nOrbitals;
	public final int[] atomicNumbers, atomOfOrb;
	public final int[][] missingOfAtom, orbsOfAtom;
	public final NDDOAtom[] atoms;
	public final NDDOOrbital[] orbitals;
	public final MoleculeInfo rm; // todo rename

	protected final SimpleMatrix B, Bforediis;
	protected final double[] earray;
	protected final transient DMatrixRMaj ddrm;

	public double energy, homo, lumo, hf, dipole;
	public double[] dipoletot;
	protected SimpleMatrix H, densityMatrix, alphaDensity, betaDensity;
	protected double ediisThreshold = config.ediis_threshold;

	protected Solution(MoleculeInfo mi, NDDOAtom[] atoms) {
		this.atoms = atoms;
		this.rm = mi;
		this.charge = mi.charge;
		this.mult = mi.mult;
		this.atomicNumbers = mi.atomicNumbers;
		this.nElectrons = mi.nElectrons;
		this.nOrbitals = mi.nOrbitals;

		earray = new double[8];
		B = new SimpleMatrix(8, 8);
		Bforediis = new SimpleMatrix(8, 8);
		ddrm = new DMatrixRMaj(nOrbitals, nOrbitals);

		orbitals = new NDDO6G[nOrbitals];
		orbsOfAtom = mi.orbsOfAtom; // variable length 2d array
		missingOfAtom = mi.missingOfAtom; // variable length 2d array, complement of orbsOfAtom
		atomOfOrb = mi.atomOfOrb; // nOrbitals length 1d array

		int overallOrbitalIndex = 0;
		for (NDDOAtom atom : atoms) {
			for (NDDOOrbital orbital : atom.getOrbitals()) {
				orbitals[overallOrbitalIndex] = orbital;
				overallOrbitalIndex++;
			}
		}
	}

	public static Solution of(MoleculeInfo mi, NDDOAtom[] atoms) {
		return of(mi, atoms, null);
	}

	public static Solution of(MoleculeInfo mi, NDDOAtom[] atoms, double[][] densityMatrices) {
		Solution s = mi.restricted ? new SolutionR(mi, atoms) : new SolutionU(mi, atoms);

		if (!mi.useEdiis) {
			s.ediisThreshold = Double.POSITIVE_INFINITY;
		}

		return s.compute(densityMatrices);
	}

	public static boolean shouldEdiis(MoleculeInfo mi, NDDOAtom[] atoms) {
		mi.getLogger().debug("Starting EDIIS test...");

		Solution stest = mi.restricted ? new SolutionR(mi, atoms) : new SolutionU(mi, atoms);

		double[] results = stest.testEdiis();
		double v = results[0] - results[1];
		boolean useEdiis = v < config.ediis_max_diff || v != v; // if Inf - Inf then EDIIS
		if (!useEdiis)
			mi.getLogger().warn("EDIIS gives higher energy ({} - {} = {})- using purely DIIS from now on.",
					results[0], results[1], v);

		mi.getLogger().debug("useEdiis = {}", useEdiis);

		return useEdiis;
	}

	public static SimpleMatrix[] findDensityMatrices(MoleculeInfo mi, NDDOAtom[] atoms) {
		Solution s = Solution.of(mi, atoms);

		if (mi.restricted) return new SimpleMatrix[]{s.densityMatrix};
		else return new SimpleMatrix[]{s.alphaDensity, s.betaDensity};
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
		List<Double> array = new LinkedList<>();

		for (double i : original.getDDRM().data) {
			array.add(i);
		}

		for (int index : tbr) {
			array.add(index, 0.0);
		}

		return new SimpleMatrix(original.numRows() + tbr.length, 1, true, Utils.toDoubles(array));
	}

	public abstract Solution withNewAtoms(NDDOAtom[] newAtoms);

	public final Solution compute() {
		compute(null);

		return this;
	}

	public final Solution compute(double[][] densityMatrices) {
		precomp();

		if (densityMatrices != null) {
			if (rm.restricted) {
				densityMatrix = new SimpleMatrix(rm.nOrbitals, rm.nOrbitals, true, densityMatrices[0]);
			}
			else {
				alphaDensity = new SimpleMatrix(rm.nOrbitals, rm.nOrbitals, true, densityMatrices[0]);
				betaDensity = new SimpleMatrix(rm.nOrbitals, rm.nOrbitals, true, densityMatrices[1]);
			}
		}

		computePrivate();

		findMatrices();
		findHf();
		findDipole();
		findHomo();

		rm.getLogger().debug("hf = {}, dipole = {}, homo = {}", hf, dipole, homo);

		return this;
	}

	private void precomp() {
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

	protected abstract void computePrivate();

	protected final SimpleMatrix ediis(int len) {
		int ediisSize = len + 1;
		int len1 = len - 1;

		SimpleMatrix mat = new SimpleMatrix(ediisSize, ediisSize);

		for (int i = 0; i < len; i++) {
			for (int j = i; j < len; j++) {
				double v = Bforediis.get(i, j);
				mat.set(i, j, v);
				mat.set(j, i, v);
			}
		}

		double[] row = new double[ediisSize];
		double[] col = new double[ediisSize];
		Arrays.fill(row, 1);
		Arrays.fill(col, 1);
		mat.setColumn(len, 0, row);
		mat.setRow(len, 0, col);
		mat.set(len, len, 0);

		SimpleMatrix rhs = SimpleMatrix.ones(ediisSize, 1);
		for (int i = 0; i < len; i++) {
			rhs.set(i, earray[i]);
		}

		double bestE = 0;
		SimpleMatrix bestEdiis = null;
		for (int i = 0; i <= len1; i++) {
			for (int[] tbr : TBRS[i]) {
				try {
					SimpleMatrix newmat = removeElements(mat, tbr);
					SimpleMatrix newrhs = removeElements(rhs, tbr);
					SimpleMatrix tempEdiis = addRows(newmat.solve(newrhs), tbr);
					tempEdiis.set(len, 0);
					boolean nonNegative = !(elementMin(tempEdiis.getDDRM()) < 0);

					if (nonNegative) {
						double e = 0;

						for (int a = 0; a < len; a++) {
							e -= earray[a] * tempEdiis.get(a);

							for (int b = 0; b < len; b++) {
								e += 0.5 * tempEdiis.get(a) * tempEdiis.get(b) * Bforediis.get(a, b);
							}
						}

						if (e < bestE) {
							bestE = e;
							bestEdiis = tempEdiis;
						}
					}
				} catch (SingularMatrixException ignored) {
				}
			}
		}

		return bestEdiis;
	}

	protected final SimpleMatrix diis(int len) {
		int ediisSize = len + 1;

		SimpleMatrix mat = new SimpleMatrix(ediisSize, ediisSize);

		for (int i = 0; i < len; i++) {
			for (int j = i; j < len; j++) {
				double v = B.get(i, j);
				mat.set(i, j, v);
				mat.set(j, i, v);
			}
		}

		double[] a = new double[ediisSize];
		Arrays.fill(a, 1);
		mat.setColumn(len, 0, a);
		mat.setRow(len, 0, a);
		mat.set(len, len, 0);

		SimpleMatrix rhs = new SimpleMatrix(ediisSize, 1);
		rhs.set(len, 0, 1);

		return mat.solve(rhs);
	}

	protected abstract void findMatrices();

	protected abstract void findHf();

	protected void findDipole() {
		double[] populations = new double[atoms.length];

		for (int j = 0; j < atoms.length; j++) {
			double sum = 0;
			for (int k : orbsOfAtom[j]) {
				sum += densityMatrix().get(k, k);
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

		double[] chargedip = new double[]{0, 0, 0};

		for (int j = 0; j < atoms.length; j++) {
			double v1 = Constants.DIPOLECONV * populations[j];
			chargedip[0] += v1 * (atoms[j].getCoordinates()[0] - com[0]);
			chargedip[1] += v1 * (atoms[j].getCoordinates()[1] - com[1]);
			chargedip[2] += v1 * (atoms[j].getCoordinates()[2] - com[2]);
		}

		double[] hybridip = new double[]{0, 0, 0};

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

	protected abstract void findHomo();

	private double[] testEdiis() {
		precomp();

		double ediishf;
		try {
			densityMatrix = alphaDensity = betaDensity = null;
			ediisThreshold = config.ediis_threshold;
			computePrivate();
			findHf();
			ediishf = hf;
		} catch (IllegalStateException e) {
			ediishf = Double.POSITIVE_INFINITY;
		}

		double noEdiishf;
		try {
			densityMatrix = alphaDensity = betaDensity = null;
			ediisThreshold = Double.POSITIVE_INFINITY;
			computePrivate();
			findHf();
			noEdiishf = hf;
		} catch (IllegalStateException e) {
			noEdiishf = Double.POSITIVE_INFINITY;
		}

		return new double[]{ediishf, noEdiishf};
	}

	public final MoleculeInfo getRm() {
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
}
