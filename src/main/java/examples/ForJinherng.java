package examples;

import frontend.TxtIO;
import nddo.NDDOAtom;
import nddo.geometry.GeometryOptimization;
import nddo.solution.Solution;
import runcycle.State;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;

import java.io.IOException;

public class ForJinherng {
	public static void main(String[] args) throws IOException {
		TxtIO.txtToText();

		RunInput input = TxtIO.readInput();

		for (RunnableMolecule molecule : input.molecules) {
			NDDOAtom[] atoms = State.getConverter().convert(molecule.atoms, input.info.npMap);

			Solution s = Solution.of(molecule, atoms);

			molecule.getLogger().info("Initial Hf: {}", s.hf);

			GeometryOptimization go = GeometryOptimization.of(s).compute();

			Solution optS = go.getS();

			molecule.getLogger().info("Final Hf: {}", optS.hf);

			molecule.getLogger().info("Dipole: {}", optS.dipole);

			molecule.getLogger().info("IE: {}", -optS.homo);
		}

		System.out.println("Press 'Enter' key to exit.");
		System.console().readLine();
	}
}
