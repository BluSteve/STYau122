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

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

public class Utils {
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
				coords[tau] = coords[tau] + 1E-7;
				perturbed[i] = atoms[i].withNewCoords(coords);
			}
			else {
				perturbed[i] = atoms[i].copy();
			}
		}

		return perturbed;
	}

	public static String getResource(String resourceName) throws IOException {
		InputStream inputStream =
				Objects.requireNonNull(Utils.class.getClassLoader().getResourceAsStream(resourceName));

		ByteArrayOutputStream result = new ByteArrayOutputStream();
		byte[] buffer = new byte[1024];
		for (int length; (length = inputStream.read(buffer)) != -1; ) {
			result.write(buffer, 0, length);
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
		int i1 = vector.numRows();
		for (int i = 0; i < i1; i++) {
			double v = vector.get(i);
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

	public static double pow(double d, double n) {
		final double absn = abs(n);
		if (absn > 16) return Math.pow(d, n);
		else if (d == 0) return 0;
		else if (d == 1) return 1;

		final double iabsn = (int) absn;
		double r = 0;

		if (absn == iabsn) { // if integral
			if (iabsn == 0) return 1;
			else if (iabsn == 1) r = d;
			else if (iabsn == 2) r = d * d;
			else if (iabsn == 3) r = d * d * d;
			else if (iabsn == 4) {
				double d2 = d * d;
				r = d2 * d2;
			}
			else if (iabsn == 5) {
				double d2 = d * d;
				r = d2 * d2 * d;
			}
			else if (iabsn == 6) {
				double d2 = d * d;
				r = d2 * d2 * d2;
			}
			else if (iabsn == 7) {
				double d2 = d * d;
				r = d2 * d2 * d2 * d;
			}
			else if (iabsn == 8) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4;
			}
			else if (iabsn == 9) {
				double d3 = d * d * d;
				r = d3 * d3 * d3;
			}
			else if (iabsn == 10) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d2;
			}
			else if (iabsn == 11) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d2 * d;
			}
			else if (iabsn == 12) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d4;
			}
			else if (iabsn == 13) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d4 * d;
			}
			else if (iabsn == 14) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d4 * d2;
			}
			else if (iabsn == 15) {
				double d2 = d * d;
				double d5 = d2 * d2 * d;
				r = d5 * d5 * d5;
			}
			else if (iabsn == 16) {
				double d2 = d * d;
				double d4 = d2 * d2;
				double d8 = d4 * d4;
				r = d8 * d8;
			}
		}
		else if (iabsn < 5){ // integer comparison faster
			if (absn == 0.5) r = sqrt(d); // hardcoded to prevent function call overhead
			else if (absn == 1.5) r = d * sqrt(d);
			else if (absn == 2.5) r = d * d * sqrt(d);
			else if (absn == 3.5) r = d * d * d * sqrt(d);
			else if (absn == 4.5) {
				double d2 = d * d;
				r = d2 * d2 * sqrt(d);
			}

			else if (absn == 0.25) r = sqrt(sqrt(d));
			else if (absn == 1.25) r = d * sqrt(sqrt(d));
			else if (absn == 2.25) r = d * d * sqrt(sqrt(d));
			else if (absn == 3.25) r = d * d * d * sqrt(sqrt(d));
			else if (absn == 4.25) {
				double d2 = d * d;
				r = d2 * d2 * sqrt(sqrt(d));
			}

			else if (absn == 0.75) {
				double root = sqrt(d);
				r = root * sqrt(root);
			}
			else if (absn == 1.75) {
				double root = sqrt(d);
				r = d * root * sqrt(root);
			}
			else if (absn == 2.75) {
				double root = sqrt(d);
				r = d * d * root * sqrt(root);
			}
			else if (absn == 3.75) {
				double root = sqrt(d);
				r = d * d * d * root * sqrt(root);
			}
			else if (absn == 4.75) {
				double d2 = d * d;
				double root = sqrt(d);
				r = d2 * d2 * root * sqrt(root);
			}

			// any higher than 4 Math.pow becomes faster.
		}
		else return Math.pow(d, n);

		if (r != 0) {
			if (n > 0) return r;
			else return 1 / r;
		}

		return Math.pow(d, n);
	}

	public static void main(String[] args) {
		double i = 1234;
		double n = -4.5;

		Random r = new Random(1);
		for (int j = 0; j < 1000; j++) {
			i = r.nextDouble();
			double pow = Utils.pow(i, n);
			double x = Math.pow(i, n) - pow;
			System.out.println("pow = " + pow);
			if (x > 1E-15) System.out.println(x);
		}
	}
}

class Pair<F extends Comparable<F>, S> implements Comparable<Pair<F, S>> {
	public final F first;
	public final S second;

	public Pair(F first, S second) {
		this.first = first;
		this.second = second;
	}

	@Override
	public int compareTo(Pair<F, S> o) {
		return this.first.compareTo(o.first);
	}

	@Override
	public String toString() {
		return "Pair{" +
				"first=" + first.toString() +
				", second=" + second.toString() +
				'}';
	}
}

