package core;

import frontend.FrontendConfig;
import nddo.solution.Solution;

import java.io.IOException;

public class Testing {
	public static void main(String[] args) throws IOException {
		FrontendConfig.init();
		Benchmarks.State state = new Benchmarks.State();
		state.setup();

		System.out.println(Solution.USE_EDIIS);
	}
}
