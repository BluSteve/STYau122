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
	private static final double coresinv = 1.0 / Runtime.getRuntime().availableProcessors();
	private static final int minBatchSize = 1;

	public static void forloop(int length, Consumer<int[]> consumer) {
		int[] inputArr = new int[length];
		for (int i = 0; i < length; i++) inputArr[i] = i;

		if (length <= minBatchSize) consumer.accept(inputArr);
		else {
			final int batchSize = Math.max((int) Math.ceil(length * coresinv), minBatchSize);
			final int batchSize2 = batchSize << 1;

			final List<RecursiveAction> subtasks = new ArrayList<>(length);

			int endSize = 0;
			for (int es = 0; endSize < length; es += batchSize) {
				int elapsedSize = es;
				endSize = length - elapsedSize < batchSize2 ? length : elapsedSize + batchSize;

				int finalEndSize = endSize;
				subtasks.add(new RecursiveAction() {
					@Override
					protected void compute() {
						consumer.accept(Arrays.copyOfRange(inputArr, elapsedSize, finalEndSize));
					}
				});
			}

			ForkJoinTask.invokeAll(subtasks);
		}
	}

	public static <T> void consume(T[] inputArr, Consumer<T[]> consumer) {
		final int length = inputArr.length;

		if (length <= minBatchSize) consumer.accept(inputArr);
		else {
			final int batchSize = Math.max((int) Math.ceil(length * coresinv), minBatchSize);
			final int batchSize2 = batchSize << 1;

			final List<RecursiveAction> subtasks = new ArrayList<>(length);

			int endSize = 0;
			for (int es = 0; endSize < length; es += batchSize) {
				int elapsedSize = es;
				endSize = length - elapsedSize < batchSize2 ? length : elapsedSize + batchSize;

				int finalEndSize = endSize;
				subtasks.add(new RecursiveAction() {
					@Override
					protected void compute() {
						consumer.accept(Arrays.copyOfRange(inputArr, elapsedSize, finalEndSize));
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
			final int batchSize = Math.max((int) Math.ceil(length * coresinv), minBatchSize);
			final int batchSize2 = batchSize << 1;

			final var results = (R[]) Array.newInstance(inputArr[0].getClass(), length);
			final List<RecursiveAction> subtasks = new ArrayList<>(length);

			int endSize = 0;
			for (int es = 0; endSize < length; es += batchSize) {
				int elapsedSize = es;
				endSize = length - elapsedSize < batchSize2 ? length : elapsedSize + batchSize;

				int finalEndSize = endSize;
				subtasks.add(new RecursiveAction() {
					@Override
					protected void compute() {
						System.arraycopy(function.apply(Arrays.copyOfRange(inputArr, elapsedSize, finalEndSize)),
								0, results, elapsedSize, finalEndSize - elapsedSize);
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
			final int batchSize = Math.max((int) Math.ceil(length * coresinv), minBatchSize);
			final int batchSize2 = batchSize << 1;

			final var results = new double[length];
			final List<RecursiveAction> subtasks = new ArrayList<>(length);

			int endSize = 0;
			for (int es = 0; endSize < length; es += batchSize) {
				int elapsedSize = es;
				endSize = length - elapsedSize < batchSize2 ? length : elapsedSize + batchSize;

				int finalEndSize = endSize;
				subtasks.add(new RecursiveAction() {
					@Override
					protected void compute() {
						System.arraycopy(function.apply(Arrays.copyOfRange(inputArr, elapsedSize, finalEndSize)),
								0, results, elapsedSize, finalEndSize - elapsedSize);
					}
				});
			}

			ForkJoinTask.invokeAll(subtasks);

			return results;
		}
	}
}
