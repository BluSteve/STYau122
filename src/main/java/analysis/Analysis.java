package analysis;

import frontend.JsonIO;
import runcycle.IMoleculeResult;
import runcycle.structs.RunOutput;

import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;

public class Analysis {
	public static void main(String[] args) throws IOException {
		RunOutput ro = JsonIO.readOutput("foranalysis");

		JsonIO.write(ro, "slimmeddown");

		IMoleculeResult[] results = ro.results.clone();
		Arrays.sort(results, Comparator.comparingDouble(IMoleculeResult::getTotalError));

		for (IMoleculeResult result : results) {
			System.out.println(result.getUpdatedRm().debugName() + ": " + result.getTotalError());
		}
	}
}
