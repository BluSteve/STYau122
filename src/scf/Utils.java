package scf;

import nddoparam.NDDOAtom;
import nddoparam.NDDOParams;
import org.apache.commons.math3.primes.Primes;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.decomposition.eig.SymmetricQRAlgorithmDecomposition_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.List;

public class Utils {
	public static final double LAMBDA = 1E-7;
	public static final double bohr = 1.88973;
	public static final int maxAtomNum = 10;

	public static double[] toDoubles(String[] strs) {
		double[] doubles = new double[strs.length];
		for (int i = 0; i < strs.length; i++) {
			doubles[i] = Double.parseDouble(strs[i]);
		}
		return doubles;
	}

	public static DoubleMatrix toDoubleMatrix(SimpleMatrix matrix) {
		double[][] array = new double[matrix.numRows()][matrix.numCols()];
		for (int r = 0; r < matrix.numRows(); r++) {
			for (int c = 0; c < matrix.numCols(); c++) {
				array[r][c] = matrix.get(r, c);
			}
		}
		return new DoubleMatrix(array);
	}

	public static boolean testEJML(DoubleMatrix x, DoubleMatrix y,
								   double limit) {
		for (int i = 0; i < y.rows; i++) {
			for (int j = 0; j < y.columns; j++) {
				if (Math.abs(Math.abs(x.get(i, j)) - Math.abs(y.get(i, j))) >
						limit) {
					System.err.println(i + ", " + j + ", " + Math.abs(
							Math.abs(x.get(i, j)) - Math.abs(y.get(i, j))));
					return false;
				}
			}
		}
		return true;
	}


	public static double[] toDoubles(List<Double> ld) {
		double[] doubles = new double[ld.size()];
		for (int i = 0; i < ld.size(); i++) {
			doubles[i] = ld.get(i);
		}
		return doubles;
	}

	public static int[] toInts(List<Integer> li) {
		int[] ints = new int[li.size()];
		for (int i = 0; i < li.size(); i++) {
			ints[i] = li.get(i);
		}
		return ints;
	}

	public static double[] bohr(double[] notbohr) {
		double[] res = new double[notbohr.length];
		for (int i = 0; i < notbohr.length; i++) res[i] = notbohr[i] * bohr;
		return res;
	}

	public static double[] debohr(double[] notbohr) {
		double[] res = new double[notbohr.length];
		for (int i = 0; i < notbohr.length; i++) res[i] = notbohr[i] / bohr;
		return res;
	}

	public static int[] findTightestTriplet(int n, int c) {
		List<Integer> primeFactors = Primes.primeFactors(n);
		primeFactors.sort((a, b) -> b - a);
		double r = Math.pow(n, 1.0 / c);
		int[] F = new int[c];
		for (int i = 0; i < c; i++) F[i] = 1;

		for (int prime : primeFactors) {
			int i;
			boolean broken = false;
			int iSmallest = 0;
			for (i = 0; i < F.length; i++) {
				int t = prime * F[i];
				if (F[i] < F[iSmallest]) iSmallest = i;
				if (t <= r) {
					F[i] = F[i] * prime;
					broken = true;
					break;
				}
			}
			if (!broken) {
				F[iSmallest] = F[iSmallest] * prime;
			}
		}
		Arrays.sort(F);
		return F;
	}

	public static int getFCores(int index) {
		int cores = Runtime.getRuntime().availableProcessors();
		int[] a = Utils.findTightestTriplet(cores, 4);
		// allocates extra cores to MoleculeRuns to account for geom opt.
		return index == 2 ? a[0] * a[3] : a[index + 1];
	}

	public static boolean hasAtomType(int[] mats, int atomType) {
		for (int p : mats) {
			if (atomType == p) {
				return true;
			}
		}
		return false;
	}

	public static double[][] to2dArray(List<List<Double>> input) {
		int size = input.size();
		double[][] array = new double[size][size];
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				array[i][j] = input.get(i).get(j);
			}
		}
		return array;
	}

	public static NDDOAtom[] perturbAtomParams(NDDOAtom[] atoms, int Z,
											   int paramNum) {
		// Z is proton number of the atom that you want to perturb the param
		// (based on paramNum) of.
		NDDOAtom[] perturbed = new NDDOAtom[atoms.length];
		try {
			Class<? extends NDDOAtom> c = atoms[0].getClass();
			Constructor<? extends NDDOAtom> ctor =
					c.getDeclaredConstructor(c,
							atoms[0].getParams().getClass());
			ctor.setAccessible(true);
			Constructor<? extends NDDOAtom> ctor2 =
					c.getDeclaredConstructor(c);
			ctor2.setAccessible(true);
			for (int i = 0; i < atoms.length; i++) {
				if (atoms[i].getAtomProperties().getZ() == Z) {
					NDDOParams params = atoms[i].getParams();
					params.modifyParam(paramNum, Utils.LAMBDA);
					perturbed[i] = ctor.newInstance(atoms[i], params);
				}
				else {
					perturbed[i] = ctor2.newInstance(atoms[i]);
				}
			}
		} catch (NoSuchMethodException | IllegalAccessException |
				InstantiationException | InvocationTargetException e) {
			e.printStackTrace();
		}
		return perturbed;
	}

	public static NDDOAtom[] perturbAtomCoords(NDDOAtom[] atoms, int atomNum,
											   int tau) {
		NDDOAtom[] perturbed = new NDDOAtom[atoms.length];
		try {
			Class<? extends NDDOAtom> c = atoms[0].getClass();
			Constructor<? extends NDDOAtom> ctor =
					c.getDeclaredConstructor(c,
							atoms[0].getCoordinates().getClass());
			ctor.setAccessible(true);
			Constructor<? extends NDDOAtom> ctor2 =
					c.getDeclaredConstructor(c);
			ctor2.setAccessible(true);
			for (int i = 0; i < atoms.length; i++) {
				if (i == atomNum) {
					double[] coords = atoms[i].getCoordinates().clone();
					coords[tau] = coords[tau] + 1E-7;
					perturbed[i] = ctor.newInstance(atoms[i], coords);
				}
				else {
					perturbed[i] = ctor2.newInstance(atoms[i]);
				}
			}
		} catch (NoSuchMethodException | IllegalAccessException |
				InstantiationException | InvocationTargetException e) {
			e.printStackTrace();
		}
		return perturbed;
	}

	public static boolean containsZ(NDDOAtom[] atoms, int Z) {
		boolean result = false;
		for (NDDOAtom atom : atoms) {
			if (atom.getAtomProperties().getZ() == Z) {
				result = true;
				break;
			}
		}
		return result;
	}

	public static int getTrainingSetSize(String trainingSet) {
		int result = 0;
		int removeH = 0;
		if (trainingSet.contains("H")) {
			result += 5;
			removeH = 1;
		}
		result += (trainingSet.length() - removeH) * 8;
		return result;
	}

	public static int numNotNull(Object[] array) {
		int count = 0;
		for (Object x : array) {
			if (x != null) count++;
		}
		return count;
	}

	public static int numIterable(int[] iterable) {

		int count = 0;

		for (int value : iterable) {

			if (value == 0) {
				count++;
			}
		}

		return count;
	}

	public static SimpleMatrix filled(int rows, int cols, int a) {
		SimpleMatrix res = new SimpleMatrix(rows, cols);
		res.fill(a);
		return res;
	}

	public static SimpleMatrix[] symEigen(SimpleMatrix sm) {
		EigenDecomposition_F64<DMatrixRMaj> evd =
				new SymmetricQRAlgorithmDecomposition_DDRM(true);
		evd.decompose(sm.copy().getDDRM());
		int noe = evd.getNumberOfEigenvalues();

		SimpleMatrix evalues = new SimpleMatrix(noe, noe);
		SimpleMatrix evectors = new SimpleMatrix(noe, noe);

		Pair<Double, DMatrixRMaj>[] epairs = new Pair[noe];

		for (int i = 0; i < noe; i++) {
			epairs[i] = new Pair<>(evd.getEigenvalue(i).real,
					evd.getEigenVector(i));
		}

		Arrays.sort(epairs);
		for (int i = 0; i < noe; i++) {
			evalues.set(i, i, epairs[i].first);
			evectors.setColumn(i, 0, epairs[i].second.data);
		}

		return new SimpleMatrix[]{evectors, evalues};
	}

	public static synchronized DoubleMatrix solve(DoubleMatrix lhs,
												  DoubleMatrix rhs) {
		return Solve.solve(lhs, rhs);
	}

	public static synchronized DoubleMatrix pinv(DoubleMatrix dm) {
		return Solve.pinv(dm);
	}

	public static SimpleMatrix[][] convertToEJML2D(
			DoubleMatrix[][] doubleMatrices) {
		SimpleMatrix[][] matrices = new SimpleMatrix[doubleMatrices.length][];
		for (int i = 0; i < doubleMatrices.length; i++) {
			DoubleMatrix[] dm = doubleMatrices[i];
			matrices[i] = new SimpleMatrix[dm.length];
			for (int j = 0; j < dm.length; j++) {
				matrices[i][j] = new SimpleMatrix(dm[j].toArray2());
			}
		}
		return matrices;
	}

	public static double[][] to2dArray(SimpleMatrix matrix) {
		double[][] array = new double[matrix.numRows()][matrix.numCols()];
		for (int r = 0; r < matrix.numRows(); r++) {
			for (int c = 0; c < matrix.numCols(); c++) {
				array[r][c] = matrix.get(r, c);
			}
		}
		return array;
	}
}

class Pair<F extends Comparable<F>, S> implements Comparable<Pair<F, S>> {
	@Override
	public int compareTo(Pair<F, S> o) {
		return this.first.compareTo(o.first);
	}

	public F first;
	public S second;

	public Pair(F first, S second) {
		this.first = first;
		this.second = second;
	}

	@Override
	public String toString() {
		return "Pair{" +
				"first=" + first.toString() +
				", second=" + second.toString() +
				'}';
	}
}
