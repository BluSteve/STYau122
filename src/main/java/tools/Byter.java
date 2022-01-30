package tools;

import java.nio.ByteBuffer;

public class Byter {
	public static byte[] toBytes(int x) {
		ByteBuffer bb = ByteBuffer.allocate(Integer.BYTES);
		bb.putInt(x);
		return bb.array();
	}

	public static byte[] toBytes(long x) {
		ByteBuffer bb = ByteBuffer.allocate(Long.BYTES);
		bb.putLong(x);
		return bb.array();
	}

	public static byte[] toBytes(float x) {
		ByteBuffer bb = ByteBuffer.allocate(Float.BYTES);
		bb.putFloat(x);
		return bb.array();
	}

	public static byte[] toBytes(double x) {
		ByteBuffer bb = ByteBuffer.allocate(Double.BYTES);
		bb.putDouble(x);
		return bb.array();
	}

	public static int toInt(byte[] b) {
		return ByteBuffer.wrap(b).getInt();
	}

	public static long toLong(byte[] b) {
		return ByteBuffer.wrap(b).getLong();

	}

	public static float toFloat(byte[] b) {
		return ByteBuffer.wrap(b).getFloat();

	}

	public static double toDouble(byte[] b) {
		return ByteBuffer.wrap(b).getDouble();
	}
}
