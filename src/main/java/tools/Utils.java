package tools;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.decomposition.eig.SymmetricQRAlgorithmDecomposition_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.simple.SimpleMatrix;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;

import static java.lang.Math.sqrt;

public class Utils {
	public static void main(String[] args) {
		double[] doubles = {37, 37, 39, 35, 40, 37, 38, 36};

		double max = 0;
		for (double v : doubles) if (v > max) max = v;
		for (int j = 0; j < doubles.length; j++) doubles[j] /= max;
		System.out.println(iqr(doubles));
		System.out.println(sd(doubles));
	}

	public static double iqr(double... values) {
		double[] valuesc = values.clone();
		Arrays.sort(valuesc);

		return valuesc[(int) (valuesc.length * 0.75)] - valuesc[(int) (valuesc.length * 0.25)];
	}

	public static double sd(double... values) {
		double mean = 0;
		for (double value : values) mean += value;
		mean /= values.length;

		double sd = 0;
		for (double value : values) {
			sd += Pow.pow(value - mean, 2);
		}
		sd /= values.length;
		sd = Math.sqrt(sd);

		return sd;
	}

	public static double[] toDoubles(String[] strs) {
		double[] doubles = new double[strs.length];
		for (int i = 0; i < strs.length; i++) {
			doubles[i] = Double.parseDouble(strs[i]);
		}
		return doubles;
	}

	public static int[] toInts(String[] strs) {
		int[] ints = new int[strs.length];
		for (int i = 0; i < strs.length; i++) {
			ints[i] = Integer.parseInt(strs[i]);
		}
		return ints;
	}

	public static double[] toDoubles(Collection<Double> ld) {
		double[] doubles = new double[ld.size()];

		int i = 0;
		for (Double aDouble : ld) {
			doubles[i] = aDouble;
			i++;
		}

		return doubles;
	}

	public static int[] toInts(Collection<Integer> li) {
		int[] ints = new int[li.size()];

		int i = 0;
		for (Integer integer : li) {
			ints[i] = integer;
			i++;
		}

		return ints;
	}

	public static boolean hasAtomType(int[] mats, int atomType) {
		for (int p : mats) {
			if (atomType == p) {
				return true;
			}
		}
		return false;
	}

	public static NDDOAtom[] perturbAtomParams(NDDOAtom[] atoms, int Z, int paramNum) {
		// Z is proton number of the atom that you want to perturb the param
		// (based on paramNum) of.
		NDDOAtom[] perturbed = new NDDOAtom[atoms.length];

		for (int i = 0; i < atoms.length; i++) {
			if (atoms[i].getAtomProperties().getZ() == Z) {
				NDDOParams params = atoms[i].getParams().copy();
				params.modifyParam(paramNum, Constants.LAMBDA);
				perturbed[i] = atoms[i].withNewParams(params);
			}
			else {
				perturbed[i] = atoms[i].copy();
			}
		}

		return perturbed;
	}

	public static NDDOAtom[] perturbAtomCoords(NDDOAtom[] atoms, int atomNum,
											   int tau) {
		NDDOAtom[] perturbed = new NDDOAtom[atoms.length];

		for (int i = 0; i < atoms.length; i++) {
			if (i == atomNum) {
				double[] coords = atoms[i].getCoordinates().clone();
				coords[tau] = coords[tau] + Constants.LAMBDA;
				perturbed[i] = atoms[i].withNewCoords(coords);
			}
			else {
				perturbed[i] = atoms[i].copy();
			}
		}

		return perturbed;
	}

	public static String getResource(String resourceName) {
		InputStream inputStream =
				Objects.requireNonNull(Utils.class.getClassLoader().getResourceAsStream(resourceName));

		ByteArrayOutputStream result = new ByteArrayOutputStream();

		try {
			byte[] buffer = new byte[1024];
			for (int length; (length = inputStream.read(buffer)) != -1; ) {
				result.write(buffer, 0, length);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		return result.toString(StandardCharsets.UTF_8);
	}

	public static List<String> getResourceAsList(String resourceName) throws IOException {
		InputStream inputStream =
				Objects.requireNonNull(Utils.class.getClassLoader().getResourceAsStream(resourceName));

		BufferedReader r = new BufferedReader(new InputStreamReader(inputStream));
		List<String> res = new ArrayList<>();

		String line;
		while ((line = r.readLine()) != null) {
			res.add(line);
		}

		return res;
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

	public static SimpleMatrix[] symEigen(SimpleMatrix sm) {
		EigenDecomposition_F64<DMatrixRMaj> evd = new SymmetricQRAlgorithmDecomposition_DDRM(true);
		evd.decompose(sm.copy().getDDRM());
		int noe = evd.getNumberOfEigenvalues();

		SimpleMatrix evalues = new SimpleMatrix(noe, noe);
		SimpleMatrix evectors = new SimpleMatrix(noe, noe);

		Pair<Double, DMatrixRMaj>[] epairs = new Pair[noe];

		for (int i = 0; i < noe; i++) {
			epairs[i] = new Pair<>(evd.getEigenvalue(i).real, evd.getEigenVector(i));
		}

		Arrays.sort(epairs);
		for (int i = 0; i < noe; i++) {
			evalues.set(i, i, epairs[i].first);
			evectors.setColumn(i, 0, epairs[i].second.data);
		}

		return new SimpleMatrix[]{evectors, evalues};
	}

	public static double mag(SimpleMatrix vector) {
		double sum = 0;
		double[] data = vector.getDDRM().data;
		for (double v : data) {
			sum += v * v;
		}

		return sqrt(sum);
	}

	public static <T> void shuffleArray(T[] array) {
		int index;
		T temp;
		Random random = new Random(Constants.RANDOM_SEED);
		for (int i = array.length - 1; i > 0; i--) {
			index = random.nextInt(i + 1);
			temp = array[index];
			array[index] = array[i];
			array[i] = temp;
		}
	}

	public static double[] bohr(double[] notbohr) {
		double[] res = new double[notbohr.length];
		for (int i = 0; i < notbohr.length; i++) res[i] = notbohr[i] * Constants.bohr;
		return res;
	}

	public static double[] debohr(double[] bohr) {
		double[] res = new double[bohr.length];
		for (int i = 0; i < bohr.length; i++) res[i] = bohr[i] / Constants.bohr;
		return res;
	}

	public static SimpleMatrix plusTrans(SimpleMatrix sm) {
		int numRows = sm.numRows();
		int numCols = sm.numCols();

		for (int i = 0; i < numRows; i++) {
			for (int j = i; j < numCols; j++) {
				double v = sm.get(i, j) + sm.get(j, i);
				sm.set(i, j, v);
				sm.set(j, i, v);
			}
		}

		return sm;
	}

	public static boolean isSimilar(SimpleMatrix x, SimpleMatrix y, double limit) {
		for (int i = 0; i < y.numRows(); i++) {
			for (int j = 0; j < y.numCols(); j++) {
				if (Math.abs(x.get(i, j) - y.get(i, j)) > limit) {
					return false;
				}
			}
		}

		return true;
	}
}
