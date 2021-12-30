package testing;

import frontend.TxtIO;
import nddo.param.ParamSecondDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;


public class Testing {
	public static void main(String[] args) throws Exception {
		RunInput input = TxtIO.readInput();
		RunnableMolecule rm = input.molecules[0];
		RunnableMolecule rm2 = input.molecules[1];

		SolutionR s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));
		SolutionU s2 = (SolutionU) Solution.of(rm2, runcycle.State.getConverter().convert(rm2.atoms, input.info.npMap));

		System.out.close();
		for (int i = 1; i < 7; i++) {
			for (int j = i; j < 7; j++) {
				boolean b = ParamSecondDerivative.verifyEquations(s, 6, i, 6, j);
				System.err.println(i + " " + j + " " + b);
				if (!b) {
					throw new IllegalArgumentException(i + " " + j);
				}
			}
		}

		System.err.println("\n\n\nUHF\n\n\n");

		for (int i = 1; i < 7; i++) {
			boolean b = ParamSecondDerivative.verifyEquations(s2, 6, i);
			System.err.println(i + " " + b);
			if (!b) {
				throw new IllegalArgumentException(i + "");
			}
		}
	}
}
