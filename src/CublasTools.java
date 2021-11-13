import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcublas.JCublas;
import jcuda.jcusolver.JCusolverDn;
import jcuda.jcusolver.cusolverDnHandle;
import jcuda.runtime.JCuda;
import org.ejml.simple.SimpleMatrix;
import scf.Utils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

import static jcuda.jcublas.cublasFillMode.CUBLAS_FILL_MODE_UPPER;
import static jcuda.jcusolver.cusolverEigMode.CUSOLVER_EIG_MODE_VECTOR;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;

public class CublasTools {
	public static void main(String[] args) {
		int n = 10000;
		int n2 = n * n;

		cusolverDnHandle handle = new cusolverDnHandle();
		JCusolverDn.cusolverDnCreate(handle);
		Random rand = new Random();

		int warmupn = 30;
		for (int i = 0; i < 5; i++) {
			SimpleMatrix warmup =
					SimpleMatrix.random_DDRM(warmupn, warmupn, 0, 1, rand);
			Utils.symEigen(warmup);
			warmup.mult(warmup);
			getSymEigenGPU(handle, warmup);
			gpuMmul(warmup, warmup);
		}

		for (n = 10; n <= 1000; n += 5) {
			SimpleMatrix a = SimpleMatrix.random_DDRM(n, n, 0, 1, rand);
			SimpleMatrix b = SimpleMatrix.random_DDRM(n, n, 0, 1, rand);

//			for (int i = 0; i < a.numRows(); i++) {
//				for (int j = i; j < a.numCols(); j++) {
//					a.set(j, i, a.get(i, j));
//				}
//			}

//			System.out.println("a = " + a);

			NanoStopWatch sw = NanoStopWatch.sw();
			gpuMmul(a, b);
			double gpu = sw.stop();
//			System.out.println("GPU = " + gpu);

			sw.start();
			a.mult(b);
			double cpu = sw.stop();
//			System.out.println("CPU = " + cpu);

			try {
				FileWriter fw =
						new FileWriter("cpuvsgpueigennomemcpy.csv", true);
				fw.write(n + "," + cpu + "," + gpu + "\n");
				fw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public static SimpleMatrix[] getSymEigenGPU(cusolverDnHandle handle,
												SimpleMatrix a) {
		// note: the eigenvalues are sorted

		int n2 = a.getNumElements();
		int n = a.numRows();
		double[] a1d = to1d(a); // input matrix in gpu format
		double[] v = new double[n2]; // eigenvectors
		double[] w = new double[n]; // eigenvalues
		int[] info_gpul = new int[1];

		Pointer h_A = Pointer.to(a1d); // h means cpu memory
		Pointer h_V = Pointer.to(v);
		Pointer h_W = Pointer.to(w);
		Pointer h_info_gpul = Pointer.to(info_gpul);
		Pointer d_A = new Pointer(); // d means gpu memory
		Pointer d_V = new Pointer();
		Pointer d_W = new Pointer();

		Pointer d_work = new Pointer();
		Pointer d_dev_info = new Pointer();

		JCuda.cudaMalloc(d_A, (long) n2 * Sizeof.DOUBLE);
		JCuda.cudaMalloc(d_W, (long) n * Sizeof.DOUBLE);
		JCuda.cudaMalloc(d_dev_info, Sizeof.INT);

		int jobz = CUSOLVER_EIG_MODE_VECTOR;
		int uplo = CUBLAS_FILL_MODE_UPPER;

		JCuda.cudaMemcpy(d_A, h_A, (long) n2 * Sizeof.DOUBLE,
				cudaMemcpyHostToDevice); // getting input to gpu

		// retrieves working space requirements
		int[] lworkl = new int[1];
		JCusolverDn.cusolverDnDsyevd_bufferSize(handle, jobz, uplo, n, d_A, n,
				d_W, lworkl);
		int lwork = lworkl[0];
		JCuda.cudaMalloc(d_work, (long) lwork * Sizeof.DOUBLE);

		JCusolverDn.cusolverDnDsyevd(handle, jobz, uplo, n, d_A, n,
				d_W, d_work, n2, d_dev_info);

		JCuda.cudaMemcpy(h_W, d_W, (long) Sizeof.DOUBLE * n,
				cudaMemcpyDeviceToHost);
		JCuda.cudaMemcpy(h_V, d_A, (long) Sizeof.DOUBLE * n2,
				cudaMemcpyDeviceToHost);
		JCuda.cudaMemcpy(h_info_gpul, d_dev_info, Sizeof.INT,
				cudaMemcpyDeviceToHost);

		JCuda.cudaFree(d_A);
		JCuda.cudaFree(d_V);
		JCuda.cudaFree(d_W);
		JCuda.cudaFree(d_work);

		SimpleMatrix eigenvalues = new SimpleMatrix(w.length, w.length);
		for (int i = 0; i < w.length; i++) {
			eigenvalues.set(i, i, w[i]);
		}

		return new SimpleMatrix[]{from1d(v), eigenvalues};
	}

	public static double[] to1d(SimpleMatrix dm) {
		double[] res = new double[dm.numRows() * dm.numRows()];
		for (int i = 0; i < dm.numRows(); i++) {
			for (int j = 0; j < dm.numRows(); j++) {
				res[i * dm.numRows() + j] = dm.get(j, i);
			}
		}
		return res;
	}

	public static SimpleMatrix from1d(double[] darray) {
		int n = (int) Math.sqrt(darray.length);
		SimpleMatrix res = new SimpleMatrix(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				res.set(j, i, darray[i * n + j]);
			}
		}
		return res;
	}

	public static double[] gpuMmul(double[] a, double[] b) {
		int n2 = a.length;
		int N = (int) Math.sqrt(n2);
		double[] c = new double[n2];

		Pointer gpuPointerA = new Pointer();
		Pointer gpuPointerB = new Pointer();
		Pointer gpuPointerC = new Pointer();

		JCublas.cublasAlloc(n2, Sizeof.DOUBLE, gpuPointerA);
		JCublas.cublasAlloc(n2, Sizeof.DOUBLE, gpuPointerB);
		JCublas.cublasAlloc(n2, Sizeof.DOUBLE, gpuPointerC);

		JCublas.cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(a), 1,
				gpuPointerA, 1);
		JCublas.cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(b), 1,
				gpuPointerB, 1);
		JCublas.cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(c), 1,
				gpuPointerC, 1);

		JCublas.cublasDgemm('n', 'n', N, N, N, 1.0, gpuPointerA, N,
				gpuPointerB, N, 0.0, gpuPointerC, N);

		JCublas.cublasGetVector(n2, Sizeof.DOUBLE, gpuPointerC, 1,
				Pointer.to(c), 1);

		JCublas.cublasFree(gpuPointerA);
		JCublas.cublasFree(gpuPointerB);
		JCublas.cublasFree(gpuPointerC);

		return c;
	}

	public static SimpleMatrix gpuMmul(SimpleMatrix a, SimpleMatrix b) {
		return from1d(gpuMmul(to1d(a), to1d(b)));
	}

}
