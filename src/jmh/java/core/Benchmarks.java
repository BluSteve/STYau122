package core;

import frontend.TxtIO;
import nddo.param.ParamHessianNew;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import org.ejml.simple.SimpleMatrix;
import org.openjdk.jmh.annotations.*;
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
	@Warmup(iterations = 3, time = 10)
	@Measurement(iterations = 3, time = 30)
	@BenchmarkMode(Mode.SampleTime)
	@OutputTimeUnit(TimeUnit.MILLISECONDS)
	public static void init(State state) {
		ParamHessianNew pg = new ParamHessianNew(state.s, state.rm.datum, null);
	}

	@org.openjdk.jmh.annotations.State(Scope.Benchmark)
	public static class State {
		public SolutionR s;
		public SimpleMatrix[] fockderivstatic;
		public SimpleMatrix x;
		public Random r = new Random(123);
		public RunnableMolecule rm;

		@Setup(Level.Trial)
		public void setup() throws IOException {
			RunInput input = TxtIO.readInput();
			rm = input.molecules[0];

			s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));

		}
	}
}
