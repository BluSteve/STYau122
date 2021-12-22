package core;

import nddo.solution.SolutionR;
import org.ejml.simple.SimpleMatrix;
import org.openjdk.jmh.annotations.*;
import tools.Pow;

import java.io.IOException;
import java.util.Random;
import java.util.concurrent.TimeUnit;

public class Benchmarks {
	public static void main(String[] args) throws Exception {
		org.openjdk.jmh.Main.main(args);
	}

	@Benchmark
	@Fork(value = 1, warmups = 0)
	@Warmup(iterations = 5, time = 5)
	@Measurement(iterations = 5, time = 5)
	@BenchmarkMode(Mode.Throughput)
	@OutputTimeUnit(TimeUnit.NANOSECONDS)
	public static void init(State state) {
		double x = 10000 * state.r.nextDouble();
		Pow.exp(x);
	}

	@org.openjdk.jmh.annotations.State(Scope.Benchmark)
	public static class State {
		public SolutionR s;
		public SimpleMatrix[] fockderivstatic;
		public SimpleMatrix x;
		public Random r = new Random(123);

		@Setup(Level.Trial)
		public void setup() throws IOException {
//			RunInput input = TxtIO.readInput();
//			RunnableMolecule rm = input.molecules[0];
//
//			s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));
//			SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s);
//			fockderivstatic = matrices[1];
//
//			System.out.println(fockderivstatic.length);
//			SimpleMatrix[] xarray = GeometrySecondDerivative.getxarrayPople(s, fockderivstatic);
//			x = xarray[0];
//			System.out.println(xarray[0]);
		}
	}
}
