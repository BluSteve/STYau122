package tools;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.function.Consumer;
import java.util.function.Function;

public final class Batcher {
	private static final double cores = Runtime.getRuntime().availableProcessors();
	private static final int minBatchSize = 1;

	public static void forloop(int length, Consumer<int[]> consumer) {
		int[] inputArr = new int[length];
		for (int i = 0; i < length; i++) inputArr[i] = i;

		if (length <= minBatchSize) consumer.accept(inputArr);
		else {
			final int batchSize = Math.max((int) Math.ceil(length / cores), minBatchSize);

			final List<RecursiveAction> subtasks = new ArrayList<>(length);

			for (int es = 0; es < length; es += batchSize) {
				int elapsedSize = es;

				subtasks.add(new RecursiveAction() {
					@Override
					protected void compute() {
						consumer.accept(Arrays.copyOfRange(inputArr, elapsedSize,
								Math.min(length, elapsedSize + batchSize)));
					}
				});
			}

			ForkJoinTask.invokeAll(subtasks);
		}
	}

	public static <T> void consume(T[] inputArr, Consumer<T[]> consumer) {
		final int length = inputArr.length;
		System.out.println("length = " + length);

		if (length <= minBatchSize) consumer.accept(inputArr);
		else {
			final int batchSize = Math.max((int) Math.ceil(length / cores), minBatchSize);

			final List<RecursiveAction> subtasks = new ArrayList<>(length);

			for (int es = 0; es < length; es += batchSize) {
				int elapsedSize = es;

				subtasks.add(new RecursiveAction() {
					@Override
					protected void compute() {
						consumer.accept(Arrays.copyOfRange(inputArr, elapsedSize,
								Math.min(length, elapsedSize + batchSize)));
					}
				});
			}

			ForkJoinTask.invokeAll(subtasks);
		}
	}

	public static <T, R> R[] apply(T[] inputArr, Function<T[], R[]> function) {
		final int length = inputArr.length;
		if (length <= minBatchSize) return function.apply(inputArr);
		else {
			final int batchSize = Math.max((int) Math.ceil(length / cores), minBatchSize);

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
	}

	public static <T> double[] applyDouble(T[] inputArr, Function<T[], double[]> function) {
		final int length = inputArr.length;
		if (length <= minBatchSize) return function.apply(inputArr);
		else {
			final int batchSize = Math.max((int) Math.ceil(length / cores), minBatchSize);

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
	}
}
