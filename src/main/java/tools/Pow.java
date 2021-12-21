package tools;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

public class Pow {
	public static final int max = 1_000_000;
	public static final double[] arr = generatearr();

	private static double[] generatearr() {
		double[] res = new double[2 * max + 1];
		int count = 0;
		for (int j = -max; j <= max; j++, count++) {
			res[count] = Math.exp(j);
		}
		return res;
	}

	public static double exp9(double x) {
		return (39916800 + x * (39916800 + x * (19958400 + x * (6652800 +
				x * (1663200 + x * (332640 + x * (55440 + x * (7920 + x * (990 + x * (110 + x * (11 + x))))))))))) *
				2.505210838544172E-8;
	}

	public static double exp(double x) {
		int xi = (int) Math.round(x);
		double xf = x - xi;
		double exi = arr[xi + max];
		if (xf == 0) return exi;
		return exi * exp9(xf);
	}

	/**
	 * Only use this instead of Math.pow() if:<br>
	 *   a. n is an integer with abs(n) <= 16<br>
	 *   b. n ends with 0.25 or 0.5 or 0.75 and abs(n) <= 4.5
	 * @param d base
	 * @param n exponent
	 * @return d^n
	 */
	public static double pow(double d, double n) {
		final double absn = abs(n);
		if (absn > 16) {
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
			else if (absn == 4.5) {
				double d2 = d * d;
				r = d2 * d2 * sqrt(d);
			}
			else if (absn == 0.75) {
				double root = sqrt(d);
				r = root * sqrt(root);
			}

			// for now unused
			else if (absn == 3.5) r = d * d * d * sqrt(d);
			else if (absn == 1.75) {
				double root = sqrt(d);
				r = d * root * sqrt(root);
			}
			else if (absn == 2.75) {
				double root = sqrt(d);
				r = d * d * root * sqrt(root);
			}
			else if (absn == 0.25) r = sqrt(sqrt(d));
			else if (absn == 1.25) r = d * sqrt(sqrt(d));
			else if (absn == 2.25) r = d * d * sqrt(sqrt(d));
			else if (absn == 3.25) r = d * d * d * sqrt(sqrt(d));
			else if (absn == 4.25) {
				double d2 = d * d;
				r = d2 * d2 * sqrt(sqrt(d));
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
			return Math.pow(d, n);
		}

		if (r != 0) {
			if (n < 0) r = 1 / r;
			return r;
		}

		return Math.pow(d, n);
	}
}
