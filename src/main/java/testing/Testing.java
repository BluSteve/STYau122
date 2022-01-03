package testing;

import frontend.TxtIO;
import nddo.param.*;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;


public class Testing {
	public static void main(String[] args) throws Exception {
		RunInput input = TxtIO.readInput();
		RunnableMolecule rm = input.molecules[0];
		RunnableMolecule rm2 = input.molecules[1];

		SolutionR s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));
//		SolutionR se = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.expGeom,
//				input.info.npMap));
		SolutionU s2 = (SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.atoms,
				input.info.npMap));
//		SolutionU se2 =
//				(SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.expGeom, input.info.npMap));

		ParamHessianNew pg = new ParamHessianNew(s, rm.datum, null);
		ParamHessianNew pg2 = new ParamHessianNew(s2, rm2.datum, null);

//		verify(s, rm.datum, null);
//		verify(s2, rm2.datum, null);

//		verifyEquations(s, 6, 6);
//		verifyEquations(s2, 6, 6);
	}

	private static void verifyEquations(Solution s, int Z1, int Z2) {
		for (int i = 1; i < 7; i++) {
			for (int j = i; j < 7; j++) {
				boolean b2 = s instanceof SolutionR ?
						ParamSecondDerivative.verifyEquations((SolutionR) s, Z1, i, Z2, j)
						: ParamSecondDerivative.verifyEquations((SolutionU) s, Z1, i, Z2, j);

				System.err.println(i + " " + j + " " + b2);
				if (!b2) {
					throw new IllegalArgumentException(i + " " + j);
				}
			}
		}
	}

	private static void verify(Solution s, double[] datum, Solution sExp) {
		ParamGradientNew pg = new ParamGradientNew(s, datum, sExp);
		ParamGradient pg2 = ParamGradient.of(s, datum, sExp).compute();

		System.out.println("pg " + new SimpleMatrix(pg.getTotalGradients()));
		System.out.println("pg2 " + new SimpleMatrix(pg2.getTotalGradients()));

		ParamHessianNew ph = new ParamHessianNew(pg);
		ParamHessian ph2 = ParamHessian.from(pg2).compute();
		double[][] a = ph.getHessian();
		double[][] b = ph2.getHessian();

		double maxe = 0;
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {

				double e = Math.abs(a[i][j] - b[i][j]) / a[i][j];
				if (e > 0.01) System.out.println(i + " " + j + " " + e + " " + a[i][j] + " " + b[i][j]);
				if (e > maxe) maxe = e;
			}
		}

		System.out.println("max error = " + maxe);
		System.out.println("---");
	}
}
