//import org.ejml.data.DMatrixRMaj;
//import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
//import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
//import org.ejml.simple.SimpleMatrix;
//import org.jblas.DoubleMatrix;
//import org.ojalgo.matrix.Primitive64Matrix;
//import org.ojalgo.random.Weibull;
//import scf.Utils;
//
//import java.util.Random;
//
//public class EJMLTesting {
//	public static void main(String[] args) {
//		Primitive64Matrix.Factory matrixFactory = Primitive64Matrix.FACTORY;
//		int n = 100;
//		DoubleMatrix jblas = new DoubleMatrix(n,n);
//		DoubleMatrix jblas2 = new DoubleMatrix(n,n);
//		SimpleMatrix ejml = new SimpleMatrix(n,n);
//		SimpleMatrix ejml2 = new SimpleMatrix(n,n);
//		Random r =new Random(123);
//		int warmupn = 30;
//		for (int i = 0; i < 5; i++) {
//			DoubleMatrix warmup = DoubleMatrix.rand(warmupn, warmupn);
//			Utils.symEigen(warmup);
//			SimpleMatrix sm = SimpleMatrix.random_DDRM(warmupn,warmupn,0,1,new Random(1234));
//			sm.mult(sm);
//			EigenDecomposition_F64<DMatrixRMaj>
//					eig = DecompositionFactory_DDRM.eig(n*n,true,true);
//			eig.decompose(sm.getDDRM());
//
//			Primitive64Matrix a = matrixFactory.makeFilled(warmupn,warmupn,new Weibull(5,2));
//			Primitive64Matrix b = matrixFactory.makeFilled(warmupn,warmupn,new Weibull(5,2));
//			a.multiply(b);
//			Testing.gpuMmul(warmup, warmup);
//		}
//		Primitive64Matrix.DenseReceiver matrixBuilder = matrixFactory.makeDense(n,n);
//		Primitive64Matrix.DenseReceiver matrixBuilder2 = matrixFactory.makeDense(n,n);
//		for (int i = 0; i < n; i++) {
//			for (int j = i; j < n; j++) {
//				double toPut = r.nextDouble();
//				jblas.put(i*n+j,toPut);
//				jblas.put(j*n+i,toPut);
//				ejml.set(i,j,toPut);
//				ejml.set(j,i,toPut);
//				matrixBuilder.set(i,j,toPut);
//				matrixBuilder.set(j,i,toPut);
//			}
//		}
//		for (int i = 0; i < n; i++) {
//			for (int j = i; j < n; j++) {
//				double toPut = r.nextDouble();
//				jblas2.put(i*n+j,toPut);
//				jblas2.put(j*n+i,toPut);
//				ejml2.set(i,j,toPut);
//				ejml2.set(j,i,toPut);
//				matrixBuilder2.set(i,j,toPut);
//				matrixBuilder2.set(j,i,toPut);
//			}
//		}
////		pp(jblas);
////		System.out.println("ejml = " + ejml);
////		BasicLogger.debug("ojsimpson = ", matrixBuilder.get());
//		System.out.println();
//
//
//		Primitive64Matrix m1 = matrixBuilder.get();
//		Primitive64Matrix m2 = matrixBuilder2.get();
//		NanoStopWatch sw = NanoStopWatch.sw();
//
//		Primitive64Matrix m3 = m1.multiply(m2);
//		System.out.println("oj = " + sw.stop());
////				BasicLogger.debug("ojsimpson = ", m3);
//sw.start();
//Testing.gpuMmul(jblas, jblas2);
//		System.out.println("gpu = " + sw.stop());
//		sw.start();
////		DoubleMatrix[] jblaseigen = Eigen.symmetricEigenvectors(jblas);
//		SimpleMatrix sm = ejml.mult(ejml2);
//
//		System.out.println("ejml = " + sw.stop());
////		pp(dm);
//
//
////		EigenDecomposition_F64<DMatrixRMaj>
////				eig = DecompositionFactory_DDRM.eig(n*n,true,true);
////		DMatrixRMaj ddrm = ejml.getDDRM();
//		sw.start();
//		DoubleMatrix dm = jblas.mmul(jblas2);
//
////		eig.decompose(ddrm);
////		DMatrixRMaj[] eigens = getAllEigen(eig);
//		System.out.println("jblas = " + sw.stop());
////		System.out.println("sm = " + sm);;
////		System.out.println("eig = " + Arrays.toString(eigens));
//
//
//
//	}
//
//	public static void pp(DoubleMatrix dm) {
//		for (int i = 0; i < dm.rows; i++) {
//			for (int j = 0; j < dm.columns; j++) {
//				System.out.printf("%25s", dm.get(i, j));
//			}
//			System.out.println();
//		}
//		System.out.println();
//	}
//
//	public static DMatrixRMaj[] getAllEigen(EigenDecomposition_F64<DMatrixRMaj> eig) {
//		int size =  eig.getNumberOfEigenvalues();
//		DMatrixRMaj[] ddrm = new DMatrixRMaj[size];
//		for (int i = 0; i < size; i++) {
//			ddrm[i] = eig.getEigenVector(i);
//		}
//		return ddrm;
//	}
//}
