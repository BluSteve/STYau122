import org.jblas.DoubleMatrix;

import java.io.IOException;
import java.util.List;

public class SyncedDoubleMatrix extends DoubleMatrix {
	public SyncedDoubleMatrix(int newRows, int newColumns, double... newData) {
		super(newRows, newColumns, newData);
	}

	public SyncedDoubleMatrix(int newRows, int newColumns) {
		super(newRows, newColumns);
	}

	public SyncedDoubleMatrix() {
	}

	public SyncedDoubleMatrix(int len) {
		super(len);
	}

	public SyncedDoubleMatrix(double[] newData) {
		super(newData);
	}

	public SyncedDoubleMatrix(String filename) throws IOException {
		super(filename);
	}

	public SyncedDoubleMatrix(double[][] data) {
		super(data);
	}

	public SyncedDoubleMatrix(List<Double> data) {
		super(data);
	}

	@Override
	public synchronized SyncedDoubleMatrix mmul(DoubleMatrix other) {
		return staticMmul(this, other);
	}

	@Override
	public synchronized DoubleMatrix mmul(double v) {
		return staticMmul(this, v);
	}

	private static synchronized SyncedDoubleMatrix staticMmul(DoubleMatrix dm1,
															  DoubleMatrix dm2) {
		return (SyncedDoubleMatrix) dm1.mmul(dm2);
	}

	private static synchronized SyncedDoubleMatrix staticMmul(DoubleMatrix dm1,
															  double v) {
		return (SyncedDoubleMatrix) dm1.mmul(v);
	}

	public static SyncedDoubleMatrix eye(int n) {
		SyncedDoubleMatrix m = new SyncedDoubleMatrix(n, n);

		for(int i = 0; i < n; ++i) {
			m.put(i, i, 1.0D);
		}

		return m;
	}

	public static SyncedDoubleMatrix ones(int rows, int columns) {
		SyncedDoubleMatrix m = new SyncedDoubleMatrix(rows, columns);

		for(int i = 0; i < rows * columns; ++i) {
			m.put(i, 1.0D);
		}

		return m;
	}

	public static SyncedDoubleMatrix zeros(int rows, int columns) {
		return new SyncedDoubleMatrix(rows, columns);
	}
}
