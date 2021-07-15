import nddoparam.GeometryOptimizationR;
import nddoparam.Solution;
import nddoparam.SolutionR;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import nddoparam.param.ParamGradient;
import org.apache.commons.lang3.time.StopWatch;
import runcycle.MoleculeRun;
import runcycle.MoleculeRunR;
import scf.AtomHandler;
import scf.Utils;

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

//        Solution expsoln = new SolutionU(exp, 0,1 );
		Solution expsoln = new SolutionR(exp, 0);
		GeometryOptimizationR opt =
				new GeometryOptimizationR(exp1, 0);
//        GeometryOptimizationU opt = new
//        GeometryOptimizationU(exp1, 0, 1);
		StopWatch sw = new StopWatch();

		ParamGradient g;


//        g = new ParamGradientR((SolutionR) opt.s, "a", datum,
//        (SolutionR) expsoln, true);
//        g = new ParamGradientU((SolutionU) opt.s, "a",
//        datum, (SolutionU) expsoln, false);
//        g.computeGradient();
//        System.out.println("Test HF (A) Derivs: " + Arrays.deepToString(g
//        .getHFDerivs
//        ()));
//        System.out.println("Geom Derivs: " + Arrays.deepToString(g
//        .getGeomDerivs()));

//        g = new ParamGradientR(opt.s, "b", datum, expsoln, true);
//        g.computeDerivs();
//        System.out.println("Test HF (B) Derivs: " + Arrays.deepToString(g
//        .getHFDerivs
//        ()));
//        System.out.println("Test Dipole (B) Derivs: " + Arrays.deepToString(g
//        .getDipoleDerivs()));
//
//        g = new ParamGradientR(opt.s, "c", datum, expsoln, true);
//        g.computeDerivs();
//        System.out.println("Test HF (C) Derivs: " + Arrays.deepToString(g
//        .getHFDerivs
//        ()));
//        System.out.println("Test IE (C) Derivs: " + Arrays.deepToString(g
//        .getIEDerivs
//        ()));
//
////        double time = 0;
////        for (int x = 0; x < 100; x++) {
////            sw.start();
////            g = new ParamGradientR((SolutionR)opt.s, "d",
// datum, (SolutionR)expsoln, true);
////            g.computeDerivs();
//////            ParamHessianR hessian = new ParamHessianR(opt.s,
// "b", datum, expsoln);
//////            hessian.computeHessian();
////            sw.stop();
////            time += sw.getTime();
////            sw.reset();
////        }
//        sw.start();
//        g = new ParamGradientR((SolutionR) opt.s, "d", datum,
//        (SolutionR) expsoln, true);
////        g = new ParamGradientU((SolutionU) opt.s, "d",
// datum, (SolutionU) expsoln, true);
//        g.computeGradient();
//
//        sw.stop();
//        System.out.println("D Time taken: " + sw.getTime());
//        System.out.println("Test HF (D) Derivs: " + Arrays.deepToString(g
//        .getHFDerivs
//        ()));
//        System.out.println("Test Dipole (D) Derivs: " + Arrays.deepToString(g
//        .getDipoleDerivs()));
//        System.out.println("Test IE (D) Derivs: " + Arrays.deepToString(g
//        .getIEDerivs
//        ()));
//        System.out.println("Test Total (D) Derivs: " + Arrays.deepToString(g
//        .getTotalGradients()));
//        g.setAnalytical(false);
//        g.computeGradient();
//        System.out.println("Test HF (D) Derivs: " + Arrays.deepToString(g
//        .getHFDerivs
//        ()));
//        System.out.println("Test Dipole (D) Derivs: " + Arrays.deepToString(g
//        .getDipoleDerivs()));
//        System.out.println("Test IE (D) Derivs: " + Arrays.deepToString(g
//        .getIEDerivs
//        ()));
//        System.out.println("Test Total (D) Derivs: " + Arrays.deepToString(g
//        .getTotalGradients()));
//        System.out.println(g.getAnalyticalError());
//
//        sw.reset();
//        sw.start();
//        ParamHessianR hessian = new ParamHessianR(
//        (SolutionR) opt.s, "d", datum, (SolutionR) expsoln,
//        true);
////        ParamHessianR hessian2 = new ParamHessianR(
// (SolutionR) opt.s, "a", datum, null, false);
////        hessian.setAnalytical(true);
//        hessian.computeHessian();
//////        ParamHessianU hessian = new ParamHessianU(
// (SolutionU) opt.s, "b", datum, (SolutionU) expsoln);
//////        hessian.computeHessian();
//        sw.stop();
//        System.out.println("Hessian time taken: " + sw.getTime());
////        System.out.println("Test Hessian: " + Arrays.deepToString(hessian
// .getHessianUnpadded()));
//
//
//        hessian2.computeHessian();
//        System.out.println(hessian.getAnalyticalError());
//        System.out.println("Test Hessian: " + Arrays.toString(hessian
//        .getHessianUT()));
//        System.out.println("Test Hessian: " + Arrays.toString(hessian2
//        .getHessianUT()));
		MoleculeRun m = null;
		double time = 0;
		for (int x = 0; x < 1; x++) {
			sw.start();

			m = new MoleculeRunR(exp, 0, exp1, datum, true, "d",
					new int[]{1, 6});
			time += sw.getTime();
			sw.stop();
			sw.reset();
		}

		System.out.println("MOLECULERUN TIME TAKEN: " + time / 1);
	}

}
