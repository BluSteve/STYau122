package nddo.param;

import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;

import java.util.Arrays;

import static nddo.param.ParamSecondDerivative.Gderiv2static;
import static nddo.param.ParamSecondDerivative.Hderiv2;

public class ParamHessianNew implements IParamHessian {
	final ParamGradientNew pg;
	final Solution s, sExp;
	final double[] datum;
	final boolean rhf, hasDip, hasIE, hasGeom;
	final int nAtomTypes, nParams;

	final double[][] hessian;
	final SimpleMatrix[][][][][] Fstatic;

	public ParamHessianNew(ParamGradientNew pg) {
		this.pg = pg;

		s = pg.s;
		rhf = pg.rhf;
		SolutionR sr = rhf ? (SolutionR) s : null;
		SolutionU su = !rhf ? (SolutionU) s : null;

		sExp = pg.s;
		datum = pg.datum;
		hasDip = pg.hasDip;
		hasIE = pg.hasIE;
		hasGeom = pg.hasGeom;

		nAtomTypes = pg.nAtomTypes;
		nParams = pg.nParams;
		int nanp = nAtomTypes * nParams;
		int[] mats = s.rm.mats;
		int[][] mnps = s.rm.mnps;

		hessian = new double[nanp][nanp];
		Fstatic = new SimpleMatrix[nAtomTypes][nAtomTypes][nParams][nParams][];

		int[][] flat = new int[nanp * nanp][4]; // Z1, Z2, param1, param2;

		int i = 0;
		for (int ZI1 = 0; ZI1 < nAtomTypes; ZI1++) {
			for (int p1 : mnps[ZI1]) {
				for (int ZI2 = 0; ZI2 < nAtomTypes; ZI2++) {
					for (int p2 : mnps[ZI2]) {
						flat[i][0] = ZI1;
						flat[i][1] = ZI2;
						flat[i][2] = p1;
						flat[i][3] = p2;
						i++;
					}
				}
			}
		}

		flat = Arrays.copyOfRange(flat, 0, i);

		for (int[] ints : flat) {
			int ZI1 = ints[0];
			int ZI2 = ints[1];
			int Z1 = mats[ZI1];
			int Z2 = mats[ZI2];
			int p1 = ints[2];
			int p2 = ints[3];

			SimpleMatrix Hderiv2 = Hderiv2(s, Z1, p1, Z2, p2);

			if (rhf) {
				Fstatic[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
						Hderiv2.plusi(Gderiv2static(sr, Z1, p1, Z2, p2))
				};

				System.out.println(
						Arrays.toString(ints) + Fstatic[ZI1][ZI2][p1][p2][0].elementSum());
			}
			else {
				SimpleMatrix[] Gderiv2static = Gderiv2static(su, Z1, p1, Z2, p2);

				Fstatic[ZI1][ZI2][p1][p2] = new SimpleMatrix[]{
						Gderiv2static[0].plusi(Hderiv2), Gderiv2static[1].plusi(Hderiv2)
				};
			}
		}

//		System.out.println("Fstatic = " + Arrays.deepToString(Fstatic));
	}

	@Override
	public Solution getS() {
		return s;
	}

	@Override
	public double[][] getHessian() {
		return hessian;
	}
}
