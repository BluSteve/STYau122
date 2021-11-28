package tools;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.am1.AM1Params;
import nddo.mndo.MNDOParams;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.layout.PatternLayout;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.decomposition.eig.SymmetricQRAlgorithmDecomposition_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.simple.SimpleMatrix;
import runcycle.input.RawInput;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
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

	public static double[] debohr(double[] bohr) {
		double[] res = new double[bohr.length];
		for (int i = 0; i < bohr.length; i++) res[i] = bohr[i] / Utils.bohr;
		return res;
	}

	public static boolean hasAtomType(int[] mats, int atomType) {
		for (int p : mats) {
			if (atomType == p) {
				return true;
			}
		}
		return false;
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

	public static double[][] to2dArray(SimpleMatrix matrix) {
		double[][] array = new double[matrix.numRows()][matrix.numCols()];
		for (int r = 0; r < matrix.numRows(); r++) {
			for (int c = 0; c < matrix.numCols(); c++) {
				array[r][c] = matrix.get(r, c);
			}
		}
		return array;
	}

	public static String getHash(String str) {
		MessageDigest digest = null;
		try {
			digest = MessageDigest.getInstance("SHA-1");
		} catch (NoSuchAlgorithmException ignored) {
		}

		assert digest != null;
		byte[] b = digest.digest(str.getBytes(StandardCharsets.UTF_8));
		ByteBuffer wrapped = ByteBuffer.wrap(b);
		long v = wrapped.getLong();

		// gets last 41 bits of hash, 36^8 is 41.34
		// ensures low collision probability up to 10k runs
		return StringUtils.leftPad(
				Long.toUnsignedString(v & 0x1FFFFFFFFFFL, 36)
						.toUpperCase(), 8, "0");
	}

	public static double mag(SimpleMatrix vector) {
		double sum = 0;
		for (int i = 0; i < vector.numRows(); i++) {
			sum += vector.get(i) * vector.get(i);
		}

		return Math.sqrt(sum);
	}

	public static void orthogonalise(SimpleMatrix[] vectors) {
		for (int i = 0; i < vectors.length; i++) {
			for (int j = 0; j < i; j++) {
				vectors[i] = vectors[i].minus(vectors[j]
						.scale(vectors[i].dot(vectors[j]) /
								vectors[j].dot(vectors[j])));
			}
		}
	}

	// todo make this async
	public static void addFileAppender(Logger logger, String filename) {
		org.apache.logging.log4j.core.Logger l =
				(org.apache.logging.log4j.core.Logger) logger;
		LoggerContext lc = l.getContext();
		FileAppender fa = FileAppender.newBuilder()
				.setName(filename)
				.withFileName("logs/" + filename + ".log")
				.setLayout(
						PatternLayout.newBuilder().withPattern(
										"%d{ISO8601} %08r %-5level " +
												"%logger{36} - %msg%n")
								.build())
				.setConfiguration(lc.getConfiguration()).build();
		fa.start();
		lc.getConfiguration().addAppender(fa);
		((org.apache.logging.log4j.core.Logger) logger).addAppender(
				lc.getConfiguration().getAppender(fa.getName()));
		lc.updateLoggers();
	}

	public static void main(String[] args) {
		System.out.println(getHash("asdf"));
	}

	public static NDDOParams[] convertToNDDOParams(RawInput ri) {
		NDDOParams[] nddoParams = null;
		switch (ri.model) {
			case "mndo":
				nddoParams = new MNDOParams[ri.params.nddoParams.length];
				for (int i = 0; i < ri.params.nddoParams.length; i++)
					nddoParams[i] = new MNDOParams(ri.params.nddoParams[i]);
				break;
			case "am1":
				nddoParams = new AM1Params[ri.params.nddoParams.length];
				for (int i = 0; i < ri.params.nddoParams.length; i++)
					nddoParams[i] = new AM1Params(ri.params.nddoParams[i]);
				break;
		}

		assert nddoParams != null;
		return nddoParams;
	}
}

class Pair<F extends Comparable<F>, S> implements Comparable<Pair<F, S>> {
	public F first;
	public S second;

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
