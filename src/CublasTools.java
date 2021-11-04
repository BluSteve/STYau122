import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcublas.JCublas;
import jcuda.jcusolver.JCusolverDn;
import jcuda.jcusolver.cusolverDnHandle;
import jcuda.runtime.JCuda;
import org.jblas.DoubleMatrix;
import scf.Utils;

import java.io.FileWriter;
import java.io.IOException;

import static jcuda.jcublas.cublasFillMode.CUBLAS_FILL_MODE_UPPER;
import static jcuda.jcusolver.cusolverEigMode.CUSOLVER_EIG_MODE_VECTOR;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;

public class CublasTools {
	public static void main(String[] args) {
//		JCuda.setExceptionsEnabled(true);
//		JCusolver.setExceptionsEnabled(true);
		int n = 10000;
		int n2 = n * n;

		cusolverDnHandle handle = new cusolverDnHandle();
		JCusolverDn.cusolverDnCreate(handle);

		int warmupn = 30;
		for (int i = 0; i < 5; i++) {
			DoubleMatrix warmup = DoubleMatrix.rand(warmupn, warmupn);
			Utils.symEigen(warmup);
			getSymEigenGPU(handle, warmup);
		}

		for ( n = 10; n <= 500; n+=2) {
			DoubleMatrix a = DoubleMatrix.rand(n, n);
			for (int i = 0; i < a.rows; i++) {
				for (int j = i; j < a.columns; j++) {
					a.put(i * a.rows + j, a.get(i, j));
				}
			}

			NanoStopWatch sw = NanoStopWatch.sw();
			getSymEigenGPU(handle, a);
			double gpu = sw.stop();
//			System.out.println("GPU = " + gpu);

			sw.start();
			DoubleMatrix[] es = Utils.symEigen(a);
			double cpu = sw.stop();
//			System.out.println("CPU = " + cpu);

			try {
				FileWriter fw = new FileWriter("cpuvsgpueigennomemcpy.csv", true);
				fw.write(n + "," + cpu + "," + gpu + "\n");
				fw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public static DoubleMatrix getSymEigenGPU(cusolverDnHandle handle,
									  DoubleMatrix a) {
		int n2 = a.length;
		int n = a.rows;
		double[] a1d = to1d(a);
		double[] v = new double[n2];
		double[] w = new double[n];
		int[] info_gpul = new int[1];

		Pointer h_A = Pointer.to(a1d);
		Pointer h_V = Pointer.to(v);
		Pointer h_W = Pointer.to(w);
		Pointer h_info_gpul = Pointer.to(info_gpul);
		Pointer d_A = new Pointer();
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
				cudaMemcpyHostToDevice);

		int[] lworkl = new int[1];
		JCusolverDn.cusolverDnDsyevd_bufferSize(handle, jobz, uplo, n, d_A, n,
				d_W, lworkl);
		int lwork = lworkl[0];
		JCuda.cudaMalloc(d_work, (long) lwork * Sizeof.DOUBLE);

//		NanoStopWatch sw = NanoStopWatch.sw();
		JCusolverDn.cusolverDnDsyevd(handle, jobz, uplo, n, d_A, n,
				d_W, d_work, n2, d_dev_info);
//		double gpu = sw.stop();

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

		return from1d(v);
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

	public static double[] to1d(DoubleMatrix dm) {
		double[] res = new double[dm.rows * dm.rows];
		for (int i = 0; i < dm.rows; i++) {
			for (int j = 0; j < dm.rows; j++) {
				res[i * dm.rows + j] = dm.get(j, i);
			}
		}
		return res;
	}

	public static DoubleMatrix from1d(double[] darray) {
		int n = (int) Math.sqrt(darray.length);
		DoubleMatrix res = new DoubleMatrix(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				res.put(j, i, darray[i * n + j]);
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

	public static DoubleMatrix gpuMmul(DoubleMatrix a, DoubleMatrix b) {
		return from1d(gpuMmul(to1d(a), to1d(b)));
	}

}
