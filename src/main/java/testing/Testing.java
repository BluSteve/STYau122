package testing;

import frontend.TxtIO;
import nddo.param.IParamGradient;
import nddo.param.ParamGradient;
import nddo.param.ParamGradientNew;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;

import java.util.Arrays;


public class Testing {
	public static void main(String[] args) throws Exception {
		RunInput input = TxtIO.readInput();
		RunnableMolecule rm = input.molecules[0];
		RunnableMolecule rm2 = input.molecules[1];

		SolutionR s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));
		SolutionR se = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.expGeom, input.info.npMap));
		SolutionU s2 = (SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.atoms, input.info.npMap));
		SolutionU se2 = (SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.expGeom, input.info.npMap));

		IParamGradient pg = new ParamGradientNew(s, rm.datum, se);
		IParamGradient pg2 = ParamGradient.of(s, rm.datum, se).compute();

		System.out.println("pg " + Arrays.deepToString(pg.getGeomDerivs()));
		System.out.println("pg2 " + Arrays.deepToString(pg2.getGeomDerivs()));

		pg = new ParamGradientNew(s2, rm2.datum, se2);
		pg2 = ParamGradient.of(s2, rm2.datum, se2).compute();

		System.out.println("pg " + Arrays.deepToString(pg.getGeomDerivs()));
		System.out.println("pg2 " + Arrays.deepToString(pg2.getGeomDerivs()));

//		SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s);
//		SimpleMatrix[] fockderivstatic = matrices[1];
//		System.out.println(GeometrySecondDerivative.hessianRoutine(s, fockderivstatic));
//
//		matrices = GeometryDerivative.gradientRoutine(s2);
//		fockderivstatic = matrices[1];
//		SimpleMatrix[] fockderivstatic2 = matrices[2];
//		System.out.println(GeometrySecondDerivative.hessianRoutine(s2, fockderivstatic, fockderivstatic2));
//
////		System.out.close();
//		for (int i = 1; i < 7; i++) {
//			for (int j = i; j < 7; j++) {
//				boolean b = ParamSecondDerivative.verifyEquations(s, 6, i, 6, j);
//				System.err.println(i + " " + j + " " + b);
//				if (!b) {
//					throw new IllegalArgumentException(i + " " + j);
//				}
//			}
//		}
//
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
