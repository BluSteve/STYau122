package testing;

import frontend.TxtIO;
import nddo.param.ParamGradient;
import nddo.param.ParamGradientNew;
import nddo.param.ParamHessian;
import nddo.param.ParamHessianNew;
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
		SolutionR se = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.expGeom, input.info.npMap));
		SolutionU s2 = (SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.atoms, input.info.npMap));
		SolutionU se2 = (SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.expGeom, input.info.npMap));

		verify(s, rm.datum, null);

		verify(s2, rm2.datum, null);

//		for (int i = 1; i < 7; i++) {
//			for (int j = i; j < 7; j++) {
//				boolean b2 = ParamSecondDerivative.verifyEquations(s2, 6, i, 6, j);
//				System.err.println(i + " " + j + " " + b2);
//				if (!b2) {
//					throw new IllegalArgumentException(i + " " + j);
//				}
//			}
//		}

//		System.err.println("\n\n\nUHF\n\n\n");
//
//		for (int i = 1; i < 7; i++) {
//			boolean b = ParamSecondDerivative.verifyEquations(s2, 6, i);
//			System.err.println(i + " " + b);
//			if (!b) {
//				throw new IllegalArgumentException(i + "");
//			}
//		}
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

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {

				double e = Math.abs(a[i][j] - b[i][j]);
				if (e > 50) System.out.println(i + " " + j + " " + e + " " + a[i][j] + " " + b[i][j]);
			}
		}

		System.out.println("max error = " + new SimpleMatrix(a).minus(new SimpleMatrix(b)).elementMaxAbs());
		System.out.println("---");
	}
}
