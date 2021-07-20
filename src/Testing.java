import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.jcublas.JCublas;
import jcuda.runtime.JCuda;
import nddoparam.GeometryOptimization;
import nddoparam.Solution;
import nddoparam.SolutionR;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import scf.AtomHandler;
import scf.Utils;

import java.io.FileWriter;
import java.io.IOException;

public class Testing {
	public static void main(String[] args) {
		try {
			testOther();
		}catch(Exception e) {
			e.printStackTrace();
		}
	}

	private static void testOther() throws IOException {
		JCuda.setExceptionsEnabled(true);
//		int N = 10000;
		for (int i = 0; i < 5; i++) {
//				DoubleMatrix.rand(10, 10).mmul(DoubleMatrix.rand(10, 10));
			gpuMmul(DoubleMatrix.rand(5000, 5000), DoubleMatrix.rand(5000, 5000));
		}
		for (int N = 1000; N <= 10000; N+=1000) {
			System.out.println("N = " + N);
			DoubleMatrix adm = DoubleMatrix.rand(N, N);
			DoubleMatrix bdm = DoubleMatrix.rand(N, N);

			StopWatch sw = new StopWatch();
//			sw.start();
//			DoubleMatrix cdm = adm.mmul(bdm);
//			sw.stop();
//			long cpu = sw.getTime();
//			System.out.println("CPU: " + cpu);
//
//			sw.reset();
			sw.start();
			DoubleMatrix cdmgpu = gpuMmul(adm, bdm);
			sw.stop();
			long gpu = sw.getTime();
			System.out.println("GPU: " + gpu);
			System.out.println();
			FileWriter fw = new FileWriter("cpuvsgpu.csv",true);
			fw.write(N + "," + -1 + "," + gpu + "\n");
			fw.close();
		}

//		System.out.println(cdm);
//		System.out.println(cdmgpu);
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
		Solution sr = new SolutionR(exp1, 0);
		GeometryOptimization opt = GeometryOptimization.of(sr);



		/*int[] atomTypes = new int[]{1, 6};

		StopWatch sw = new StopWatch();
		RawMolecule rm = new RawMolecule();
		rm.mats = atomTypes;
		rm.mnps = new int[][] {MNDOParams.T1ParamNums, MNDOParams.T2ParamNums};
		SolutionR expsoln = (new SolutionR(exp, 0)).setRm(rm);
//		SolutionU expsoln = (new SolutionU(exp, 0,1 )).setRm(rm);

		Solution S = opt.s.setRm(rm);
		sw.start();
		ParamGradient G = ParamGradient.of(S, datum, expsoln).compute();
		ParamHessian H = ParamHessian.from(G).compute();
		sw.stop();
		System.err.println(Arrays.deepToString(H.getHessian()));
//		System.err.println(
//				Arrays.deepToString(H.getHessian(new int[]{1, 5, 6},
//						new int[][]{MNDOParams.T1ParamNums,
//								MNDOParams.T2ParamNums,
//								MNDOParams.T2ParamNums})));
		System.err.println(sw.getTime());*/
	}
}
