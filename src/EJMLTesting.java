import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import scf.Utils;

import java.util.Arrays;
import java.util.Random;

public class EJMLTesting {
	public static void main(String[] args) {
		int n = 3;
		DoubleMatrix jblas = new DoubleMatrix(n, n);
		DoubleMatrix jblas2 = new DoubleMatrix(n, n);
		SimpleMatrix ejml = new SimpleMatrix(n, n);
		SimpleMatrix ejml2 = new SimpleMatrix(n, n);
		Random r = new Random(123);
		int warmupn = 30;
		for (int i = 0; i < 5; i++) {
			DoubleMatrix warmup = DoubleMatrix.rand(warmupn, warmupn);
			SimpleMatrix sm = SimpleMatrix.random_DDRM(warmupn, warmupn, 0, 1,
					new Random(1234));
			sm.mult(sm);
			EigenDecomposition_F64<DMatrixRMaj>
					eig = DecompositionFactory_DDRM.eig(n * n, true, true);
			eig.decompose(sm.getDDRM());

		}
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				double toPut = r.nextDouble();
				jblas.put(i * n + j, toPut);
				jblas.put(j * n + i, toPut);
				ejml.set(i, j, toPut);
				ejml.set(j, i, toPut);
			}
		}
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				double toPut = r.nextDouble();
				jblas2.put(i * n + j, toPut);
				jblas2.put(j * n + i, toPut);
				ejml2.set(i, j, toPut);
				ejml2.set(j, i, toPut);
			}
		}

		SimpleMatrix sm = new SimpleMatrix(new double[][] {
				new double[]{1,3,2},
				new double[]{3,4,5},
				new double[]{2,5,6}
		});
		DoubleMatrix dm = new DoubleMatrix(new double[][] {
				new double[]{1,3,2},
				new double[]{3,4,5},
				new double[]{2,5,6}
		});

		System.out.println("ejml = " + sm);

		System.out.println(Arrays.toString(Utils.symEigen(sm)));
		System.out.println(sm.normF());
		System.out.println(dm.norm2());
		System.out.println("CommonOps_DDRM.elementMax(sm = " +
				CommonOps_DDRM.elementMax(sm.getDDRM()));

		pp(Eigen.symmetricEigenvectors(dm)[0]);
		pp(Eigen.symmetricEigenvectors(dm)[1]);


	}

	public static void pp(DoubleMatrix dm) {
		for (int i = 0; i < dm.rows; i++) {
			for (int j = 0; j < dm.columns; j++) {
				System.out.printf("%25s", dm.get(i, j));
			}
			System.out.println();
		}
		System.out.println();
	}

	public static DMatrixRMaj[] getAllEigen(
			EigenDecomposition_F64<DMatrixRMaj> eig) {
		int size = eig.getNumberOfEigenvalues();
		DMatrixRMaj[] ddrm = new DMatrixRMaj[size];
		for (int i = 0; i < size; i++) {
			ddrm[i] = eig.getEigenVector(i);
		}
		return ddrm;
	}
}
