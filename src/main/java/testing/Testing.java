package testing;

import frontend.TxtIO;
import nddo.geometry.GeometryDerivative;
import nddo.math.PopleThiel;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import org.ejml.simple.SimpleMatrix;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;


public class Testing {
	public static void main(String[] args) throws Exception {
		RunInput input = TxtIO.readInput();
		RunnableMolecule rm = input.molecules[0];

		SolutionR s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));
		SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s);
		SimpleMatrix[] fockderivstatic = matrices[1];

		SimpleMatrix[] xarray = PopleThiel.pople(s, fockderivstatic);
		System.out.println(xarray[0]);
	}
}
