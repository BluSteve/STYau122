package testing;

import frontend.FrontendConfig;
import frontend.TxtIO;
import nddo.Config;
import nddo.geometry.GeometryDerivative;
import nddo.math.PopleThiel;
import nddo.param.*;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;


public class Testing {
	public static void main(String[] args) throws Exception {
		Config config = FrontendConfig.config;
		RunInput input = TxtIO.readInput();
		RunnableMolecule rm = input.molecules[0];
		RunnableMolecule rm2 = input.molecules[1];

		SolutionR s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));
		SolutionR se = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.expGeom,
				input.info.npMap));
		SolutionU s2 = (SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.atoms,
				input.info.npMap));
		SolutionU se2 =
				(SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.expGeom, input.info.npMap));

		SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s2);
		SimpleMatrix[] falpha = matrices[1];
		SimpleMatrix[] fbeta = matrices[2];

		System.out.println(PopleThiel.verify(s, GeometryDerivative.gradientRoutine(s)[1], null, 1e-5));
		System.out.println(PopleThiel.verify(s2, falpha, fbeta, 1e-5));
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

	private static void verifyEquations(SolutionU su, int Z1) {
		for (int i = 1; i < 7; i++) {
			boolean b2 = ParamSecondDerivative.verifyEquations(su, Z1, i);
			System.err.println(i + " " + b2);
			if (!b2) {
				throw new IllegalArgumentException(i + "");
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
		int maxi = 0, maxj = 0;
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {
				double e = Math.abs(a[i][j] - b[i][j]) / a[i][j];
				if (e > 0.01) System.out.printf("%2d %2d %f %f %f%n", i, j, e, a[i][j], b[i][j]);
				if (e > maxe) {
					maxe = e;
					maxi = i;
					maxj = j;
				}
			}
		}

		System.out.println(maxi + " " + maxj + ": max error = " + maxe);
		System.out.println("---");
	}
}
