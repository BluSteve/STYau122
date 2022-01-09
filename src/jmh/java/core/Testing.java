package core;

import frontend.FrontendConfig;

import java.io.IOException;

public class Testing {
	public static void main(String[] args) throws IOException {
		FrontendConfig.init();
		Benchmarks.State state = new Benchmarks.State();
		state.setup();

		state.s.testEdiis();
	}
}
