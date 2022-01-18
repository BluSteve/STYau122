package testing;

import frontend.FrontendConfig;
import frontend.TxtIO;
import nddo.param.*;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;
import runcycle.structs.Serializer;


public class Testing {
	public static void main(String[] args) throws Exception {
		FrontendConfig.init();
		RunInput input = TxtIO.readInput("molecules.txt");
		RunIterator iterator = new RunIterator(input, 1);
		RunOutput ro = iterator.next();
		System.out.println("Serializer.gson.toJson(ro) = " + Serializer.gson.toJson(ro));

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
