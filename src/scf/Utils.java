package scf;

import nddoparam.NDDOAtom;
import nddoparam.NDDOParams;
import org.apache.commons.math3.primes.Primes;
import org.ejml.simple.SimpleEVD;
import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.jblas.Solve;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.stream.Collectors;

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

	public static DoubleMatrix toDoubleMatrix(SimpleMatrix matrix)
	{
		double[][] array = new double[matrix.numRows()][matrix.numCols()];
		for (int r=0; r<matrix.numRows(); r++)
		{
			for (int c=0; c<matrix.numCols(); c++)
			{
				array[r][c] = matrix.get(r,c);
			}
		}
		return new DoubleMatrix (array);
	}

	public static double[] vectorToDoubleArray (SimpleMatrix mat) {
		double[] arr = new double[mat.numRows()];
		for (int i = 0; i < mat.numRows(); i++) {
			arr[i] = mat.get(i, 0);
		}

		return arr;
	}



	public static SimpleMatrix[] SymmetricEigenvalueDecomposition(SimpleMatrix mat) {//massive TODO

		SimpleEVD<SimpleMatrix> eig = mat.eig();

		double[] unsortedeigenvalues = new double[mat.numRows()];

		SimpleMatrix eigenvalueMatrix = new SimpleMatrix(mat.numRows(), 1);

		for (int i = 0; i < mat.numRows(); i++) {
			if (eig.getEigenvalue(i).isReal()) {
				unsortedeigenvalues[i] = eig.getEigenvalue(i).real;
			}
			else {
				System.err.println ("eigenvalues aren't real!");
				System.exit(0);
			}
		}


		int[] indexes = new int[mat.numRows()];

		for (int i = 0; i < indexes.length; i++) {
			indexes[i] = i;
		}

		List<Integer> indexesCopy = Arrays.stream(indexes).boxed().collect(Collectors.toList());
		List<Integer> copyofindexes = new ArrayList<Integer>();

		for (int i: indexesCopy) {
			copyofindexes.add(i);
		}
		ArrayList<Integer> sortedList = new ArrayList<Integer>(indexesCopy);
		Collections.sort(sortedList, Comparator.comparing(s -> unsortedeigenvalues[copyofindexes.indexOf(s)]));

		SimpleMatrix eigenvectorMatrix = new SimpleMatrix(mat.numRows(), 0);


		for (int i = 0; i < mat.numRows(); i++) {
			eigenvalueMatrix.set(i, 0, unsortedeigenvalues[sortedList.get(i)]);
			eigenvectorMatrix = eigenvectorMatrix.concatColumns(eig.getEigenVector(sortedList.get(i)));
		}

		return new SimpleMatrix[] {eigenvectorMatrix, eigenvalueMatrix};
	}

	public static boolean testEJML (DoubleMatrix x, DoubleMatrix y,double limit) {
		for (int i = 0; i < y.rows; i++) {
			for (int j = 0; j < y.columns; j++) {
				if (Math.abs(Math.abs(x.get(i, j)) - Math.abs(y.get(i, j))) > limit) {
					System.err.println (i + ", " + j + ", " + Math.abs(Math.abs(x.get(i, j)) - Math.abs(y.get(i, j))));
					return false;
				}
			}
		}
		return true;
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

	public static int numNotNull(DoubleMatrix[] rarray) {
		int count = 0;
		for (DoubleMatrix r : rarray) {
			if (r != null) {
				count++;
			}
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

	public static synchronized DoubleMatrix[] symEigen(DoubleMatrix dm) {
		return Eigen.symmetricEigenvectors(dm);
	}
	
	public static synchronized DoubleMatrix solve(DoubleMatrix lhs, DoubleMatrix rhs) {
		return Solve.solve(lhs, rhs);
	}

	public static synchronized DoubleMatrix pinv(DoubleMatrix dm) {
		return Solve.pinv(dm);
	}
}
