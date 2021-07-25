public class NanoStopWatch {
	long time;
	static final double DIVISOR = 1E6;

	private NanoStopWatch() {
		start();
	}

	public static NanoStopWatch sw() {
		return new NanoStopWatch();
	}

	public double stop() {
		return (System.nanoTime() - time) / DIVISOR;
	}

	public double start() {
		time = System.nanoTime();
		return time / DIVISOR;
	}
}
