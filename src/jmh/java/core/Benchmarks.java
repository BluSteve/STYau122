package core;

import nddo.NDDOParams;
import nddo.geometry.GeometryDerivative;
import nddo.geometry.GeometrySecondDerivative;
import nddo.param.ParamDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import org.ejml.simple.SimpleMatrix;
import org.openjdk.jmh.annotations.*;
import runcycle.input.InputHandler;
import runcycle.input.RawInput;
import runcycle.input.RawMolecule;
import tools.Utils;

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
			RawInput ri = InputHandler.processInput("subset");
			RawMolecule rm = ri.molecules[0];

			NDDOParams[] npMap = Utils.getNpMap(ri);
			s = (SolutionR) Solution.of(rm, Utils.toNDDOAtoms(ri.model, rm.atoms, npMap));
			SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s);
			fockderivstatic = matrices[1];

			System.out.println(fockderivstatic.length);
			SimpleMatrix[] x = ParamDerivative.xArrayLimitedPople(s, fockderivstatic);
			System.out.println(x[0]);
		}
	}
}
