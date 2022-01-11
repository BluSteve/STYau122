package core;

import frontend.TxtIO;
import nddo.solution.Solution;
import org.openjdk.jmh.annotations.*;
import runcycle.structs.InputInfo;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;

import java.io.IOException;
import java.util.Random;
import java.util.concurrent.TimeUnit;

public class Benchmarks {
	public static void main(String[] args) throws Exception {
		org.openjdk.jmh.Main.main(args);
	}

	@Benchmark
	@Fork(value = 1, warmups = 0)
	@Warmup(iterations = 3, time = 5)
	@Measurement(iterations = 3, time = 5)
	@BenchmarkMode(Mode.SampleTime)
	@OutputTimeUnit(TimeUnit.MICROSECONDS)
	public static void init(State state) {
		Solution s = Solution.of(state.rm, runcycle.State.getConverter().convert(state.rm.atoms, state.info.npMap));
	}

	@org.openjdk.jmh.annotations.State(Scope.Benchmark)
	public static class State {
		public Random r = new Random(123);
		public Solution s;
		public RunnableMolecule rm;
		public InputInfo info;

		@Setup(Level.Trial)
		public void setup() throws IOException {
			RunInput input = TxtIO.readInput("molecules.txt");
			RunnableMolecule rm = input.molecules[1];
			System.out.println("rm = " + rm);
			this.rm = rm;
			this.info = input.info;

			s = Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));

			System.out.println("s.hf = " + s.hf);
			System.out.println("s.dipole = " + s.dipole);
			System.out.println("s.homo = " + s.homo);
		}
	}
}
