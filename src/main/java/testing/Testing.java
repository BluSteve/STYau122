package testing;

import frontend.TxtIO;
import nddo.param.ParamSecondDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;


public class Testing {
	public static void main(String[] args) throws Exception {
		RunInput input = TxtIO.readInput();
		RunnableMolecule rm = input.molecules[0];

		SolutionR s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));

		System.err.println(ParamSecondDerivative.verifyEquations(s, 6, 1, 6, 1));
	}
}
