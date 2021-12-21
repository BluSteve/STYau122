package tools;

import java.util.Random;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

public class Pow {
	private static final int max = 250; // anything outside is 0 or infinity
	private static final int interp = 4;
	private static final double interpinv = 1.0 / interp;
	private static final int total = interp * max;
	private static final double[] arr = generatearr();

	public static void main(String[] args) {
		double i = interpinv / 2;
		System.out.println("worst case = " + Math.abs(Math.exp(i) - Pow.exp(i)));

		Random r = new Random(123);
		double error = 0;
		for (int j = 0; j < 100000; j++) {
			i = r.nextGaussian();
			error += Math.abs(Math.exp(i) - Pow.exp(i));
		}
		System.out.println("  avg case = " + error / 100000);
	}

	private static double[] generatearr() {
		double[] res = new double[2 * interp * max + 1];
		int count = 0;
		for (int j = -total; j <= total; j++, count++) {
			res[count] = Math.exp(j * interpinv);
		}
		return res;
	}

	public static double exp8(double x) {
		return (3628800 + x * (3628800 + x * (1814400 +
				x * (604800 + x * (151200 + x * (30240 + x * (5040 + x * (720 + x * (90 + x * (10 + x)))))))))) *
				2.755731922398589E-7;
	}

	public static double exp(double x) { // 0.25
		if (x < -max) return 0;
		if (x > max) return Double.POSITIVE_INFINITY;
		int xi = (int) Math.round(x * interp); // 1
		double xf = x - xi * interpinv; // 0.25
		double exi = arr[xi + total]; // 41
		if (xf == 0) return exi;
		return exi * exp8(xf);
	}

	/**
	 * Only use this instead of Math.pow() if:<br>
	 * a. n is an integer with abs(n) <= 16<br>
	 * b. n ends with 0.25 or 0.5 or 0.75 and abs(n) <= 4.5
	 *
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
