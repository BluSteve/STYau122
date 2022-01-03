package testing;

import frontend.TxtIO;
import nddo.param.*;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import org.ejml.simple.SimpleMatrix;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;

import java.util.Arrays;


public class Testing {
	public static void main(String[] args) throws Exception {
		RunInput input = TxtIO.readInput();
		RunnableMolecule rm = input.molecules[0];
		RunnableMolecule rm2 = input.molecules[1];

		SolutionR s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));
//		SolutionR se = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.expGeom, input.info.npMap));
//		SolutionU s2 = (SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.atoms, input.info.npMap));
//		SolutionU se2 = (SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.expGeom, input.info.npMap));

		ParamGradientNew pg = new ParamGradientNew(s, rm.datum, null);
		ParamGradient pg2 = ParamGradient.of(s, rm.datum, null).compute();

		System.out.println("pg " + Arrays.deepToString(pg.getIEDerivs()));
		System.out.println("pg2 " + Arrays.deepToString(pg2.getIEDerivs()));
		ParamHessianNew ph = new ParamHessianNew(pg);
		double[][] a = ph.getHessian();

		ParamHessian ph2 = ParamHessian.from(pg2).compute();
		double[][] b = ph2.getHessian();

		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a[0].length; j++) {
				System.out.println(i + " " + j + " " + Math.abs(a[i][j] - b[i][j]) + " " + a[i][j] + " " + b[i][j]);
			}
		}

		System.out.println(new SimpleMatrix(a).minus(new SimpleMatrix(b)).elementMaxAbs());

//		pg = new ParamGradientNew(s2, rm2.datum, se2);
////		pg2 = ParamGradient.of(s2, rm2.datum, se2).compute();
//
//		System.out.println("pg " + Arrays.deepToString(pg.getHfDerivs()));
//		System.out.println("pg2 " + Arrays.deepToString(pg2.getTotalGradients()));
//		ph = new ParamHessianNew(pg);
//
//		SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s);
//		SimpleMatrix[] fockderivstatic = matrices[1];
//		System.out.println(GeometrySecondDerivative.hessianRoutine(s, fockderivstatic));
//
//		matrices = GeometryDerivative.gradientRoutine(s2);
//		fockderivstatic = matrices[1];
//		SimpleMatrix[] fockderivstatic2 = matrices[2];
//		System.out.println(GeometrySecondDerivative.hessianRoutine(s2, fockderivstatic, fockderivstatic2));

//		System.out.close();
		for (int i = 1; i < 7; i++) {
			for (int j = i; j < 7; j++) {
				boolean b2 = ParamSecondDerivative.verifyEquations(s, 7, i, 7, j);
				System.err.println(i + " " + j + " " + b2);
				if (!b2) {
					throw new IllegalArgumentException(i + " " + j);
				}
			}
		}

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
}
