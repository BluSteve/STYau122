import nddoparam.GeometryOptimizationR;
import nddoparam.Solution;
import nddoparam.SolutionR;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamGradient;
import nddoparam.param.ParamHessian;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.input.RawMolecule;
import scf.AtomHandler;
import scf.Utils;

import java.util.Arrays;

public class Testing {
	public static void main(String[] args) {
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
		System.out.close();
		GeometryOptimizationR opt = new GeometryOptimizationR(exp1, 0);
//        GeometryOptimizationU opt = new GeometryOptimizationU(exp1, 0, 1);
		int[] atomTypes = new int[]{1, 6};

		StopWatch sw = new StopWatch();
		sw.start();
		RawMolecule rm = new RawMolecule();
		rm.mats = atomTypes;
		rm.mnps = new int[][] {MNDOParams.T1ParamNums, MNDOParams.T2ParamNums};
		SolutionR expsoln = (new SolutionR(exp, 0)).setRm(rm);
//		SolutionU expsoln = (new SolutionU(exp, 0,1 )).setRm(rm);

		Solution S = opt.s.setRm(rm);
		ParamGradient G = ParamGradient.of(S, datum, expsoln).compute();
		ParamHessian H = ParamHessian.from(G).compute();
		sw.stop();
		System.err.println(Arrays.deepToString(H.getHessianRaw()));
//		System.err.println(
//				Arrays.deepToString(H.getHessian(new int[]{1, 5, 6},
//						new int[][]{MNDOParams.T1ParamNums,
//								MNDOParams.T2ParamNums,
//								MNDOParams.T2ParamNums})));
		System.err.println(sw.getTime());
	}
}
