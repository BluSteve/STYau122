package tools;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

public class Pow {
	private static final Map<Double, Double> cache = new ConcurrentHashMap<>(1048576);
	public static long hits = 0;
	public static long total = 0;

	private static double szudzik(double a, double b) {
		return a >= b ? a * a + a + b : a + b * b;
	}

	public static double pow(double d, double n) {
		total++;
		double szudzik = szudzik(d, n);
		if (cache.containsKey(szudzik)) {
			hits++;
			return cache.get(szudzik);
		}

		final double absn = abs(n);
		if (absn > 16) {
			double r = Math.pow(d, n);
			cache.put(szudzik, r);
			return Math.pow(d, n);
		}
		else if (d == 0) return 0;
		else if (d == 1) return 1;

		final double iabsn = (int) absn;
		double r = 0;

		if (absn == iabsn) { // if integral
			if (iabsn == 0) return 1;
			else if (iabsn == 1) r = d;
			else if (iabsn == 2) r = d * d;
			else if (iabsn == 3) r = d * d * d;
			else if (iabsn == 4) {
				double d2 = d * d;
				r = d2 * d2;
			}
			else if (iabsn == 5) {
				double d2 = d * d;
				r = d2 * d2 * d;
			}
			else if (iabsn == 6) {
				double d2 = d * d;
				r = d2 * d2 * d2;
			}
			else if (iabsn == 7) {
				double d2 = d * d;
				r = d2 * d2 * d2 * d;
			}
			else if (iabsn == 8) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4;
			}
			else if (iabsn == 9) {
				double d3 = d * d * d;
				r = d3 * d3 * d3;
			}
			else if (iabsn == 10) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d2;
			}
			else if (iabsn == 11) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d2 * d;
			}
			else if (iabsn == 12) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d4;
			}
			else if (iabsn == 13) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d4 * d;
			}
			else if (iabsn == 14) {
				double d2 = d * d;
				double d4 = d2 * d2;
				r = d4 * d4 * d4 * d2;
			}
			else if (iabsn == 15) {
				double d2 = d * d;
				double d5 = d2 * d2 * d;
				r = d5 * d5 * d5;
			}
			else if (iabsn == 16) {
				double d2 = d * d;
				double d4 = d2 * d2;
				double d8 = d4 * d4;
				r = d8 * d8;
			}
		}
		else if (iabsn < 5) { // integer comparison faster
			if (absn == 0.5) r = sqrt(d); // hardcoded to prevent function call overhead
			else if (absn == 1.5) r = d * sqrt(d);
			else if (absn == 2.5) r = d * d * sqrt(d);
			else if (absn == 3.5) r = d * d * d * sqrt(d);
			else if (absn == 4.5) {
				double d2 = d * d;
				r = d2 * d2 * sqrt(d);
			}

			else if (absn == 0.25) r = sqrt(sqrt(d));
			else if (absn == 1.25) r = d * sqrt(sqrt(d));
			else if (absn == 2.25) r = d * d * sqrt(sqrt(d));
			else if (absn == 3.25) r = d * d * d * sqrt(sqrt(d));
			else if (absn == 4.25) {
				double d2 = d * d;
				r = d2 * d2 * sqrt(sqrt(d));
			}

			else if (absn == 0.75) {
				double root = sqrt(d);
				r = root * sqrt(root);
			}
			else if (absn == 1.75) {
				double root = sqrt(d);
				r = d * root * sqrt(root);
			}
			else if (absn == 2.75) {
				double root = sqrt(d);
				r = d * d * root * sqrt(root);
			}
			else if (absn == 3.75) {
				double root = sqrt(d);
				r = d * d * d * root * sqrt(root);
			}
			else if (absn == 4.75) {
				double d2 = d * d;
				double root = sqrt(d);
				r = d2 * d2 * root * sqrt(root);
			}

			// any higher than 4 Math.pow becomes faster.
		}
		else {
			r = Math.pow(d, n);
			cache.put(szudzik, r);
			return r;
		}

		if (r != 0) {
			if (n < 0) r = 1 / r;
			cache.put(szudzik, r);
			return r;
		}

		r = Math.pow(d, n);
		cache.put(szudzik, r);
		return r;
	}
}
