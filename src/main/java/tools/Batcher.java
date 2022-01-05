package tools;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;

public final class Batcher {
	private static final double coresinv = 1.0 / Runtime.getRuntime().availableProcessors();
	private static final int minBatchSize = 1;
	private static final int minBatchSize2 = minBatchSize << 1;

	public static <T> void consume(T[] inputArr, Consumer<T[]> consumer) {
		final int length = inputArr.length;

		if (length <= minBatchSize2) consumer.accept(inputArr);
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

	public static <T> T[] apply(T[] inputArr, Function<T[], T[]> function) { // assumes output can cast into input
		final T[] outputArr = (T[]) Array.newInstance(inputArr[0].getClass(), inputArr.length);

		apply(inputArr, outputArr, function);

		return outputArr;
	}

	public static <T, R> R[] apply(T[] inputArr, Class<R[]> rClass, Function<T[], R[]> function) {
		final R[] outputArr = (R[]) Array.newInstance(rClass.getComponentType(), inputArr.length);

		apply(inputArr, outputArr, function);

		return outputArr;
	}

	private static <T, R> void apply(T[] inputArr, R[] outputArr, Function<T[], R[]> function) {
		final int length = inputArr.length;
		if (length <= minBatchSize2) System.arraycopy(function.apply(inputArr), 0, outputArr, 0, length);
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
						System.arraycopy(function.apply(Arrays.copyOfRange(inputArr, elapsedSize, finalEndSize)),
								0, outputArr, elapsedSize, finalEndSize - elapsedSize);
					}
				});
			}

			ForkJoinTask.invokeAll(subtasks);
		}
	}

	public static <T, U, R> R[] apply(T[] inputArr, U[] inputArr2, Class<R[]> rClass,
									   BiFunction<T[], U[], R[]> function) {
		final R[] outputArr = (R[]) Array.newInstance(rClass.getComponentType(), inputArr.length);

		apply(inputArr, inputArr2, outputArr, function);

		return outputArr;
	}

	private static <T, U, R> void apply(T[] inputArr, U[] inputArr2, R[] outputArr,
										BiFunction<T[], U[], R[]> function) {
		final int length = inputArr.length;
		if (length <= minBatchSize2) System.arraycopy(function.apply(inputArr, inputArr2), 0, outputArr, 0, length);
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
						System.arraycopy(
								function.apply(
										Arrays.copyOfRange(inputArr, elapsedSize, finalEndSize),
										Arrays.copyOfRange(inputArr2, elapsedSize, finalEndSize)
								), 0, outputArr, elapsedSize, finalEndSize - elapsedSize);
					}
				});
			}

			ForkJoinTask.invokeAll(subtasks);
		}
	}

	public static <T> double[] applyDouble(T[] inputArr, Function<T[], double[]> function) {
		final int length = inputArr.length;
		if (length <= minBatchSize2) return function.apply(inputArr);
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
