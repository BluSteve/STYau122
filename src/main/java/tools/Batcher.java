package tools;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.function.Function;

public class Batcher {
	private static final double cores = Runtime.getRuntime().availableProcessors();

	public static <T, R> R[] apply(T[] inputArr, Function<T[], R[]> function) {
		final int length = inputArr.length;
		final int batchSize = Math.max((int) Math.ceil(length / cores), 1);

		final var results = (R[]) Array.newInstance(inputArr[0].getClass(), length);
		final List<RecursiveAction> subtasks = new ArrayList<>(length);

		for (int es = 0; es < length; es += batchSize) {
			int elapsedSize = es;

			subtasks.add(new RecursiveAction() {
				@Override
				protected void compute() {
					int endSize = Math.min(length, elapsedSize + batchSize);

					System.arraycopy(function.apply(Arrays.copyOfRange(inputArr, elapsedSize, endSize)),
							0, results, elapsedSize, endSize - elapsedSize);
				}
			});
		}

		ForkJoinTask.invokeAll(subtasks);

		return results;
	}

	public static <T> double[] applyDouble(T[] inputArr, Function<T[], double[]> function) {
		final int length = inputArr.length;
		final int batchSize = Math.max((int) Math.ceil(length / cores), 1);
		
		final var results = new double[length];
		final List<RecursiveAction> subtasks = new ArrayList<>(length);

		for (int es = 0; es < length; es += batchSize) {
			int elapsedSize = es;

			subtasks.add(new RecursiveAction() {
				@Override
				protected void compute() {
					int endSize = Math.min(length, elapsedSize + batchSize);

					System.arraycopy(function.apply(Arrays.copyOfRange(inputArr, elapsedSize, endSize)),
							0, results, elapsedSize, endSize - elapsedSize);
				}
			});
		}

		ForkJoinTask.invokeAll(subtasks);

		return results;
	}

	public static <T> int[] applyInt(T[] inputArr, Function<T[], int[]> function) {
		final int length = inputArr.length;
		final int batchSize = Math.max((int) Math.ceil(length / cores), 1);

		final var results = new int[length];
		final List<RecursiveAction> subtasks = new ArrayList<>(length);

		for (int es = 0; es < length; es += batchSize) {
			int elapsedSize = es;

			subtasks.add(new RecursiveAction() {
				@Override
				protected void compute() {
					int endSize = Math.min(length, elapsedSize + batchSize);

					System.arraycopy(function.apply(Arrays.copyOfRange(inputArr, elapsedSize, endSize)),
							0, results, elapsedSize, endSize - elapsedSize);
				}
			});
		}

		ForkJoinTask.invokeAll(subtasks);

		return results;
	}

	public static <T> boolean[] applyBoolean(T[] inputArr, Function<T[], boolean[]> function) {
		final int length = inputArr.length;
		final int batchSize = Math.max((int) Math.ceil(length / cores), 1);

		final var results = new boolean[length];
		final List<RecursiveAction> subtasks = new ArrayList<>(length);

		for (int es = 0; es < length; es += batchSize) {
			int elapsedSize = es;

			subtasks.add(new RecursiveAction() {
				@Override
				protected void compute() {
					int endSize = Math.min(length, elapsedSize + batchSize);

					System.arraycopy(function.apply(Arrays.copyOfRange(inputArr, elapsedSize, endSize)),
							0, results, elapsedSize, endSize - elapsedSize);
				}
			});
		}

		ForkJoinTask.invokeAll(subtasks);

		return results;
	}
}
