package core;

import frontend.TxtIO;
import nddo.geometry.GeometryDerivative;
import nddo.geometry.GeometrySecondDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import org.ejml.simple.SimpleMatrix;
import org.openjdk.jmh.annotations.*;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;

import java.io.IOException;
import java.util.concurrent.TimeUnit;

public class Benchmarks {
	public static void main(String[] args) throws Exception {
		org.openjdk.jmh.Main.main(args);
	}

	@Benchmark
	@Fork(value = 2, warmups = 0)
	@Warmup(iterations = 5, time = 5)
	@Measurement(iterations = 5, time = 5)
	@BenchmarkMode(Mode.SampleTime)
	@OutputTimeUnit(TimeUnit.NANOSECONDS)
	public static void init(State state) {
		GeometrySecondDerivative.getxarrayPople(state.s, state.fockderivstatic);
	}

	@org.openjdk.jmh.annotations.State(Scope.Benchmark)
	public static class State {
		public SolutionR s;
		public SimpleMatrix[] fockderivstatic;

		@Setup(Level.Trial)
		public void setup() throws IOException {
			RunInput input = TxtIO.readInput();
			RunnableMolecule rm = input.molecules[0];

			s = (SolutionR) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));
			SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s);
			fockderivstatic = matrices[1];

			System.out.println(fockderivstatic.length);
			SimpleMatrix[] xarray = GeometrySecondDerivative.getxarrayPople(s, fockderivstatic);
			System.out.println(xarray[0]);
		}
	}
}
