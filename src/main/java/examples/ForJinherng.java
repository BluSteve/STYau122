package examples;

import frontend.TxtIO;
import nddo.NDDOAtom;
import nddo.geometry.GeometryOptimization;
import nddo.solution.Solution;
import runcycle.State;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;
import tools.Batcher;

import java.io.IOException;

public class ForJinherng {
	public static void main(String[] args) throws IOException {
		TxtIO.txtToText();

		RunInput input = TxtIO.readInput();

		Batcher.consume(input.molecules, subset -> {
			for (RunnableMolecule molecule : subset) {
				NDDOAtom[] atoms = State.getConverter().convert(molecule.atoms, input.info.npMap);

				Solution s = Solution.of(molecule, atoms);

				GeometryOptimization go = GeometryOptimization.of(s).compute();

				Solution optS = go.getS();

				molecule.getLogger().info("Initial Hf: {}, final Hf: {}, dipole: {}, IE: {}",
						s.hf, optS.hf, optS.dipole, -optS.homo);
			}
		});

		if (System.console() != null) {
			System.out.println("Press 'Enter' key to exit.");
			System.console().readLine();
		}
	}
}
