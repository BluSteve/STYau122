package tools;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import org.apache.commons.lang3.StringUtils;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.decomposition.eig.SymmetricQRAlgorithmDecomposition_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.simple.SimpleMatrix;

import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Arrays;
import java.util.Collection;

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
				NDDOParams params = atoms[i].getParams();
				params.modifyParam(paramNum, Constants.LAMBDA);
				perturbed[i] = atoms[i].withNewParams(params);
			}
			else {
				perturbed[i] = atoms[i].clone();
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
				perturbed[i] = atoms[i].clone();
			}
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
		EigenDecomposition_F64<DMatrixRMaj> evd = new SymmetricQRAlgorithmDecomposition_DDRM(true);
		evd.decompose(sm.copy().getDDRM());
		int noe = evd.getNumberOfEigenvalues();

		SimpleMatrix evalues = new SimpleMatrix(noe, noe);
		SimpleMatrix evectors = new SimpleMatrix(noe, noe);

		@SuppressWarnings("unchecked")
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

//	public static NDDOParams[] getNpMap(RawInput ri) { // todo move this somewhere else
//		NDDOParams[] npMap = new NDDOParams[Constants.maxAtomNum];
//
//		for (int i = 0; i < ri.atomTypes.length; i++) {
//			switch (ri.model) {
//				case MNDO:
//					npMap[ri.atomTypes[i]] = new NDDOParams(ri.params.nddoParams[i]);
//					break;
//				case AM1:
//					npMap[ri.atomTypes[i]] = new AM1Params(ri.params.nddoParams[i]);
//					break;
//			}
//		}
//
//		return npMap;
//	}
//
//	// todo you know what to do
//	public static NDDOAtom[] toNDDOAtoms(Model model, RawAtom[] ras, NDDOParams[] npMap) {
//		NDDOAtom[] atoms = new NDDOAtom[ras.length];
//
//		for (int i = 0; i < ras.length; i++) {
//			RawAtom ra = ras[i];
//			switch (model) {
//				case MNDO:
//					atoms[i] = new MNDOAtom(AtomProperties.getAtoms()[ra.Z], ra.coords, npMap[ra.Z]);
//					break;
//				case AM1:
//					atoms[i] = new AM1Atom(AtomProperties.getAtoms()[ra.Z], ra.coords,
//							(AM1Params) npMap[ra.Z]);
//					break;
//			}
//		}
//
//		return atoms;
//	}
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

