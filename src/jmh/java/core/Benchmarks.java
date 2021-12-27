package core;

import frontend.TxtIO;
import nddo.geometry.GeometryDerivative;
import nddo.geometry.GeometrySecondDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionU;
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
	@Warmup(iterations = 2, time = 5)
	@Measurement(iterations = 3, time = 5)
	@BenchmarkMode(Mode.SampleTime)
	@OutputTimeUnit(TimeUnit.MILLISECONDS)
	public static void thiel(State state) {
//		PopleThiel.thiel(state.s, state.fockderivstatic);
	}

	@org.openjdk.jmh.annotations.State(Scope.Benchmark)
	public static class State {
		public SolutionU s;
		public SimpleMatrix[] fockderivstatic, fderivalpha, fderivbeta;
		public SimpleMatrix x;
		public Random r = new Random(123);

		@Setup(Level.Trial)
		public void setup() throws IOException {
			RunInput input = TxtIO.readInput();
			RunnableMolecule rm = input.molecules[0];

			s = (SolutionU) Solution.of(rm, runcycle.State.getConverter().convert(rm.atoms, input.info.npMap));
			SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s);
			fderivalpha = matrices[1];
			fderivbeta = matrices[2];

			System.out.println(fderivalpha.length);
			System.out.println(fderivbeta.length);
			System.out.println(GeometrySecondDerivative.hessianRoutine(s, fderivalpha, fderivbeta));
//			x = xarray[0];
//			System.out.println(xarray[0]);
		}
	}
}
