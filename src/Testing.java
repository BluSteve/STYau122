import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcublas.JCublas;
import jcuda.runtime.JCuda;
import nddoparam.Solution;
import nddoparam.SolutionNew;
import nddoparam.SolutionR;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import scf.AtomHandler;
import scf.Utils;

import java.io.IOException;
import java.util.Random;

public class Testing {
	public static void main(String[] args) {
		try {
			testMain();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void testTransferSpeed() throws Exception {
		JCuda.setExceptionsEnabled(true);
		for (int i = 0; i < 5; i++) {
			gpuMmul(DoubleMatrix.rand(1000, 1000),
					DoubleMatrix.rand(1000, 1000));
		}
		Random r = new Random(123);
		int n2 = 10;
		int x = 20;
		double[][] arrays = new double[x][n2];
		double[] array = new double[x * n2];
		for (int j = 0; j < x; j++) {
			for (int i = 0; i < n2; i++) {
				arrays[j][i] = r.nextDouble();
			}
		}

		long start= System.nanoTime();
		Pointer[] pointers = new Pointer[x];
		for (int i = 0; i < pointers.length; i++){
			pointers[i] =new Pointer();
			JCublas.cublasAlloc(n2, Sizeof.DOUBLE, pointers[i]);
			JCublas.cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(arrays[i]), 1,
					pointers[i], 1);
		}
		long separate = System.nanoTime()-start;
		System.out.println("separate = " + separate/1e6);

		start = System.nanoTime();
		for (int j = 0; j < x; j++) {
			System.arraycopy(arrays[j], 0, array, j * n2, n2);
		}
		Pointer gpuPointerB = new Pointer();
		JCublas.cublasAlloc(n2*x, Sizeof.DOUBLE, gpuPointerB);
		JCublas.cublasSetVector(n2*x, Sizeof.DOUBLE, Pointer.to(array), 1,
				gpuPointerB, 1);
		pointers = new Pointer[x];
		for (int i = 0; i < pointers.length; i++){
			pointers[i] =new Pointer();
//			JCuda.cudaMemcpy(
//			JCublas.cublasZcopy(
		}
		long together = System.nanoTime()-start;
		System.out.println("together = " + together/1e6);
	}

	private static void testOther() throws IOException, InterruptedException {
		JCuda.setExceptionsEnabled(true);
		StopWatch sw = new StopWatch();
		// warmup
		for (int i = 0; i < 10; i++) {
			DoubleMatrix.rand(1000, 1000).mmul(DoubleMatrix.rand(1000));
			gpuMmul(DoubleMatrix.rand(1000, 1000),
					DoubleMatrix.rand(1000, 1000));
		}
//		for (int N = 2; N <= 1021; N += 50) {
//			System.out.println("N = " + N);
//			DoubleMatrix adm = DoubleMatrix.rand(N, N);
//			DoubleMatrix bdm = DoubleMatrix.rand(N, N);
//
//			StopWatch sw = new StopWatch();
////			sw.start();
////			DoubleMatrix cdm = adm.mmul(bdm);
////			sw.stop();
////			long cpu = sw.getTime();
////			System.out.println("CPU: " + cpu);
//
//			sw.reset();
//			sw.start();
//			DoubleMatrix cdmgpu = gpuMmul(adm, bdm);
//			sw.stop();
//			long gpu = sw.getTime();
//			System.out.println("GPU: " + gpu);
//			System.out.println();
//			FileWriter fw = new FileWriter("cpuvsgpu2.csv", true);
//			fw.write(N + "," + cpu + "," + gpu + "\n");
//			fw.close();
//			TimeUnit.SECONDS.sleep(1);
//		}

		int s = 1000;
		int n = 432;
		DoubleMatrix[] dms = new DoubleMatrix[n];
		for (int i = 0; i < dms.length; i++) {
			dms[i] = DoubleMatrix.rand(s, s);
		}

		// 2 seconds just to convert from 1d lmao, s=1000,n=432
		double[][] dms1d = new double[dms.length][];
		for (int i = 0; i < dms.length; i++) {
			dms1d[i] = to1d(dms[i]);
		}
		long start = System.nanoTime();

		// 4 seconds
		double[] result = gpuMmul(dms1d);
		long gpu = System.nanoTime() - start;

		DoubleMatrix dmresgpu = from1d(result);
		System.out.println("GPU = " + gpu / 1E6 + " ");


		DoubleMatrix dmres = dms[0];
		dmres = dmres.mmul(dms[1]);

		start = System.nanoTime();
		dmres = dms[0];
		for (int i = 1; i < dms.length; i++) {
			dmres = dmres.mmul(dms[i]);
		}
		long cpu = System.nanoTime() - start;

		System.out.println("CPU = " + cpu / 1E6 + " ");
		JCublas.cublasShutdown();
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

	public static double[] gpuMmul(double[][] arrays) {
		int n2 = arrays[0].length;
		int N = (int) Math.sqrt(n2);
		double[] c = new double[n2];

		Pointer gpuPointerA = new Pointer();
		Pointer gpuPointerB = new Pointer();
		Pointer gpuPointerC = new Pointer();
		JCublas.cublasAlloc(n2, Sizeof.DOUBLE, gpuPointerA);
		JCublas.cublasAlloc(n2, Sizeof.DOUBLE, gpuPointerB);
		JCublas.cublasAlloc(n2, Sizeof.DOUBLE, gpuPointerC);
		JCublas.cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(arrays[0]), 1,
				gpuPointerA, 1);
		JCublas.cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(arrays[1]), 1,
				gpuPointerB, 1);
		JCublas.cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(c), 1,
				gpuPointerC, 1);


		JCublas.cublasDgemm('n', 'n', N, N, N, 1.0, gpuPointerA, N,
				gpuPointerB, N, 0.0, gpuPointerC, N);
		for (int i = 2; i < arrays.length; i++) {
			Pointer gpuPointerD = new Pointer();
			JCublas.cublasAlloc(n2, Sizeof.DOUBLE, gpuPointerD);
			JCublas.cublasSetVector(n2, Sizeof.DOUBLE, Pointer.to(arrays[i])
					, 1, gpuPointerD, 1);

			JCublas.cublasDgemm('n', 'n', N, N, N, 1.0, gpuPointerC, N,
					gpuPointerD, N, 0.0, gpuPointerC, N);
			JCublas.cublasFree(gpuPointerD);
		}

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


	private static void testMain() {
		AtomHandler.populateAtoms();
		MNDOParams h = new MNDOParams(2.92397599125172,
				-6.222578482830868, 0.0, -12.200235077462583, 0.0,
				1.0693232546199132, 0.0, -13.00142320543855, 12.848,
				0.0, 0.0, 0.0, 0.0);
		MNDOParams c = new MNDOParams(2.5572499654157435, -18.854021376560777,
				-8.377666892780198, -52.57072065877964, -39.05266019981942,
				1.838438013363027, 1.805140784089995, -120.60738371097112,
				12.23, 11.47,
				2.43, 11.08, 9.84);
		MNDOAtom atom1 = new MNDOAtom(AtomHandler.atomsMap.get("H"),
				new double[]{0.635 * Utils.bohr, 0.639 * Utils.bohr,
						0.635 * Utils.bohr},
				h);
		MNDOAtom atom2 = new MNDOAtom(AtomHandler.atomsMap.get("H"),
				new double[]{0.635 * Utils.bohr, -0.635 * Utils.bohr,
						-0.639 * Utils.bohr}, h);
		MNDOAtom atom3 = new MNDOAtom(AtomHandler.atomsMap.get("H"),
				new double[]{-0.639 * Utils.bohr, -0.635 * Utils.bohr,
						0.635 * Utils.bohr}, h);
		MNDOAtom atom4 = new MNDOAtom(AtomHandler.atomsMap.get("H"),
				new double[]{-0.639 * Utils.bohr, 0.639 * Utils.bohr,
						-0.639 * Utils.bohr}, h);
		MNDOAtom carbon = new MNDOAtom(AtomHandler.atomsMap.get("C"),
				new double[]{-0.0021 * Utils.bohr, 0.0021 * Utils.bohr,
						-0.0021 * Utils.bohr}, c);

		MNDOAtom[] atoms = new MNDOAtom[]{atom1,atom2,atom3,atom4,carbon};

		MNDOAtom[] exp = new MNDOAtom[]{
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{0.6304 * Utils.bohr, 0.6304 * Utils.bohr,
								0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{0.6304 * Utils.bohr, -0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{-0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr,
								0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{-0.6304 * Utils.bohr, 0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("C"),
						new double[]{0, 0, 0}, c)};
		MNDOAtom[] exp1 = new MNDOAtom[]{
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{0.6304 * Utils.bohr, 0.6304 * Utils.bohr,
								0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{0.6304 * Utils.bohr, -0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{-0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr,
								0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{-0.6304 * Utils.bohr, 0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("C"),
						new double[]{0, 0, 0}, c)};
		double[] datum = new double[]{-17.9, 0, 13.6};
		StopWatch sw = new StopWatch();
		Solution sn =new SolutionNew(atoms, 0);;
		System.out.println(" = " );

		sw.start();
		sn = new SolutionNew(atoms, 0);
//		GeometryOptimization.of(sr).compute();
		sw.stop();

		System.out.println("sr.energy = " + sn.energy);
//		System.out.println("sw.getTime() = " + sw.getTime());
		sw.reset();
		Solution sr ;

		sw.start();
		sr = new SolutionR(atoms, 0);
//		GeometryOptimization.of(sr).compute();

		sw.stop();
		System.out.println("sr.energy = " + sr.energy);
		System.out.println("sw.getTime() = " + sw.getTime());
	}
}
