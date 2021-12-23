package nddo.scf;

import nddo.structs.OrbitalProperties;
import tools.Pow;

import java.util.Arrays;

/*
this is the GTO class, and all the relevant GTO integral routines are implemented here. The GTO object represents a
GTO of the form (x-coordinates[0])^i (y-coordinates[1])^j (z-coordinates[2])^k exp(-exponent * r^2)*/
public class GTO extends Orbital {
	public final double N;
	public final double exponent; /*radial exponent*/

	public GTO(OrbitalProperties op, double[] coordinates, double exponent) {
		super(op, coordinates);
		this.exponent = exponent;

		this.N = Pow.pow(2 / Math.PI, 0.75) * Pow.pow(2, L) *
				Pow.pow(exponent, (2 * L + 3) / 4.0) /
				Math.sqrt(fact2(2 * i - 1) * fact2(2 * j - 1) * fact2(2 * k - 1));
	}

	static int fact2(int num) {
		int sum = 1;
		for (int i = num; i > 0; i -= 2) {
			sum *= i;
		}

		return sum;
	}

	public static double R(double[] P, double[] C) {
		double val = (P[0] - C[0]) * (P[0] - C[0]) + (P[1] - C[1]) * (P[1] - C[1]) + (P[2] - C[2]) * (P[2] - C[2]);

		return Math.sqrt(val);
	}

	public static double S(GTO a, GTO b) {
		return a.getN() * b.getN() * I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
				I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
				I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
	}

	public static double Sgd(GTO a, GTO b, int tau) {
		switch (tau) {
			case 0:
				return a.getN() * b.getN() *
						Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
						I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
						I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
			case 1:
				return a.getN() * b.getN() * I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
						Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
						I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
			case 2:
				return a.getN() * b.getN() * I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
						I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
						Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
		}
		return 0;
	}

	public static double Sg2d(GTO a, GTO b, int tau1, int tau2) {
		int A = Math.min(tau1, tau2);
		int B = Math.max(tau1, tau2);

		switch (A) {
			case 0:
				switch (B) {
					case 0: /*derivative wrt x and x*/
						return a.getN() * b.getN() *
								Ig2d(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
					case 1: /* derivative wrt x and y*/
						return a.getN() * b.getN() *
								Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
					case 2: /* derivative wrt x and z*/
						return a.getN() * b.getN() *
								Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

				}
			case 1:
				switch (B) {
					case 1: /* derivative wrt y and y*/
						return a.getN() * b.getN() *
								I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								Ig2d(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
					case 2: /* derivative wrt y and z*/
						return a.getN() * b.getN() *
								I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
				}
			case 2: /* derivative wrt z and z*/
				return a.getN() * b.getN() * I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
						I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
						Ig2d(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
		}
		return 0;
	}

	public static double Salphapd(GTO a, GTO b, int type) {
		switch (type) {
			case 0:
				return a.Nalphapd() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

			case 1:
				return a.getN() * b.Nalphapd()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2]);
		}

		return 0;
	}

	public static double Salphadiagp2d(GTO a, GTO b, int type) {
		switch (type) {
			case 0:
				return a.Nalphap2d() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.Nalphapd() * b.getN()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.Nalphapd() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.Nalphapd() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphadiagp2d(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphadiagp2d(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphadiagp2d(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

			case 1:
				return a.getN() * b.Nalphap2d()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.Nalphapd()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.Nalphapd()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.Nalphapd()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphadiagp2d(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphadiagp2d(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphadiagp2d(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2]);
		}

		return 0;
	}

	public static double Salphacrossp2d(GTO a, GTO b) {
		return a.Nalphapd() * b.Nalphapd()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.Nalphapd()
				* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.Nalphapd()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.Nalphapd()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.Nalphapd() * b.getN()
				* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.Nalphapd() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.Nalphapd() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

				+ a.getN() * b.getN()
				* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

				+ a.getN() * b.getN()
				* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
				* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

				+ a.getN() * b.getN()
				* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
				* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* Ialphacrossp2d(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Ialphacrossp2d(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Ialphacrossp2d(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
	}

	public static double Salphapgd(GTO a, GTO b, int type, int tau) {
		switch (tau) {
			case 0:
				switch (type) {
					case 0:
						return a.Nalphapd() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Ialphapgd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

					case 1:
						return a.getN() * b.Nalphapd()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- a.getN() * b.getN()
								* Ialphapgd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2]);
				}
			case 1:
				switch (type) {
					case 0:
						return a.Nalphapd() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapgd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

					case 1:
						return a.getN() * b.Nalphapd()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapgd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2]);
				}
			case 2:
				switch (type) {
					case 0:
						return a.Nalphapd() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapgd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

					case 1:
						return a.getN() * b.Nalphapd()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapgd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2]);
				}
		}

		return 0;
	}

	public static double Salphacrossp2gd(GTO a, GTO b, int tau) {
		switch (tau) {
			case 0:
				return a.Nalphapd() * b.Nalphapd()
						* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.Nalphapd()
						* Ialphapgd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.Nalphapd()
						* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.Nalphapd()
						* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						- a.Nalphapd() * b.getN()
						* Ialphapgd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.Nalphapd() * b.getN()
						* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.Nalphapd() * b.getN()
						* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapgd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapgd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						- a.getN() * b.getN()
						* Ialphapgd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						- a.getN() * b.getN()
						* Ialphapgd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphacrossp2gd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphacrossp2d(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphacrossp2d(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

			case 1:
				return a.Nalphapd() * b.Nalphapd()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.Nalphapd()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.Nalphapd()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapgd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.Nalphapd()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.Nalphapd() * b.getN()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						- a.Nalphapd() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapgd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.Nalphapd() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						- a.getN() * b.getN()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapgd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* Ialphapgd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapgd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						- a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapgd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphacrossp2d(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphacrossp2gd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphacrossp2d(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

			case 2:
				return a.Nalphapd() * b.Nalphapd()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.Nalphapd()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.Nalphapd()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.Nalphapd()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapgd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.Nalphapd() * b.getN()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.Nalphapd() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						- a.Nalphapd() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapgd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						- a.getN() * b.getN()
						* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapgd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						- a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapgd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphapgd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* Ialphapgd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Ialphacrossp2d(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ialphacrossp2d(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ialphacrossp2gd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
		}

		return 0;
	}

	public static double Salphadiagp2gd(GTO a, GTO b, int type, int tau) {
		switch (tau) {
			case 0:
				switch (type) {
					case 0:
						return a.Nalphap2d() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.Nalphapd() * b.getN()
								* Ialphapgd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.Nalphapd() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.Nalphapd() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Ialphapgd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Ialphapgd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])


								+ a.getN() * b.getN()
								* Ialphadiagp2gd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphadiagp2d(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphadiagp2d(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

					case 1:
						return a.getN() * b.Nalphap2d()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- 2 * a.getN() * b.Nalphapd()
								* Ialphapgd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.Nalphapd()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.Nalphapd()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

								- 2 * a.getN() * b.getN()
								* Ialphapgd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- 2 * a.getN() * b.getN()
								* Ialphapgd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])


								- a.getN() * b.getN()
								* Ialphadiagp2gd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphadiagp2d(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* Igd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphadiagp2d(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2]);
				}

			case 1:
				switch (type) {
					case 0:
						return a.Nalphap2d() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.Nalphapd() * b.getN()
								* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.Nalphapd() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapgd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.Nalphapd() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapgd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapgd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])


								+ a.getN() * b.getN()
								* Ialphadiagp2d(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphadiagp2gd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphadiagp2d(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

					case 1:
						return a.getN() * b.Nalphap2d()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.Nalphapd()
								* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- 2 * a.getN() * b.Nalphapd()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapgd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.Nalphapd()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

								- 2 * a.getN() * b.getN()
								* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* Ialphapgd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

								- 2 * a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapgd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* Ialphapd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])


								+ a.getN() * b.getN()
								* Ialphadiagp2d(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphadiagp2gd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Igd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphadiagp2d(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2]);
				}

			case 2:
				switch (type) {
					case 0:
						return a.Nalphap2d() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.Nalphapd() * b.getN()
								* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.Nalphapd() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.Nalphapd() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapgd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Ialphapd(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapgd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapgd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])


								+ a.getN() * b.getN()
								* Ialphadiagp2d(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphadiagp2d(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphadiagp2gd(a.k, b.k, a.exponent, b.exponent,
								b.coordinates[2] - a.coordinates[2]);

					case 1:
						return a.getN() * b.Nalphap2d()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.Nalphapd()
								* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ 2 * a.getN() * b.Nalphapd()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- 2 * a.getN() * b.Nalphapd()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapgd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

								+ 2 * a.getN() * b.getN()
								* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- 2 * a.getN() * b.getN()
								* Ialphapd(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphapgd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

								- 2 * a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphapd(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* Ialphapgd(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])


								+ a.getN() * b.getN()
								* Ialphadiagp2d(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								+ a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* Ialphadiagp2d(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
								* Igd(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

								- a.getN() * b.getN()
								* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
								* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
								* Ialphadiagp2gd(b.k, a.k, b.exponent, a.exponent,
								a.coordinates[2] - b.coordinates[2]);
				}
		}

		return 0;
	}

	private static double I(int l1, int l2, double a1, double a2, double R) {
		double num = Math.sqrt(Math.PI / (a1 + a2)) * Pow.exp(-a1 * a2 * R * R / (a1 + a2));
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return num;
					case 1:
						return -a1 * R / (a1 + a2) * num;
				}
				break;
			case 1:
				switch (l2) {
					case 0:
						return a2 * R / (a1 + a2) * num;
					case 1:
						return (1 / (2 * (a1 + a2)) - a1 * a2 * R * R / ((a1 + a2) * (a1 + a2))) * num;
				}
				break;
		}

		return 0;
	}

	private static double Igd(int l1, int l2, double a1, double a2, double R) {
		double num = Math.sqrt(Math.PI / (a1 + a2)) * Pow.exp(-a1 * a2 * R * R / (a1 + a2));
		double v = 2 * a1 * a2 * R * R / (a1 + a2);
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return num * 2 * a1 * a2 * R / (a1 + a2);
					case 1:
						return a1 / (a1 + a2) * num * (1 - v);
				}
				break;
			case 1:
				switch (l2) {
					case 0:
						return a2 / (a1 + a2) * num * (v - 1);
					case 1:
						return num * a1 * a2 * R / Pow.pow(a1 + a2, 2) * (3 - 2 * R * R * a1 * a2 / (a1 + a2));
				}
				break;
		}
		return 0;
	}

	private static double Ig2d(int l1, int l2, double a1, double a2, double R) {
		double num = Math.sqrt(Math.PI / (a1 + a2)) * Pow.exp(-a1 * a2 * R * R / (a1 + a2));
		double v = 2 * a1 * a2 / (a1 + a2) * R * R;
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return num * 2 * a1 * a2 / (a1 + a2) * (v - 1);
					case 1:
						return 2 * a1 * a1 * a2 * R / ((a1 + a2) * (a1 + a2)) * num * (3 - v);
				}
				break;
			case 1:
				switch (l2) {
					case 0:
						return 2 * a1 * a2 * a2 * R / ((a1 + a2) * (a1 + a2)) * num * (v - 3);
					case 1:
						return a1 * a2 * Pow.pow(a1 + a2, -2) * num *
								(2 * R * a1 * a2 / (a1 + a2) * (3 * R - 2 * R * R * R * a1 * a2 / (a1 + a2)) +
										6 * R * R * a1 * a2 / (a1 + a2) - 3);
				}
				break;
		}
		return 0;
	}

	private static double Ialphapd(int l1, int l2, double a1, double a2, double R) {
		double exp = Pow.exp(-a1 * a2 * R * R / (a1 + a2));
		double sqrt = Math.sqrt(Math.PI / (a1 + a2));
		double derivnum = -exp * R * R *
				(a2 / (a1 + a2) - a1 * a2 / ((a1 + a2) * (a1 + a2))) * sqrt
				- exp * 0.5 * Math.sqrt(Math.PI) / Pow.pow(a1 + a2, 1.5);

		double num = sqrt * exp;

		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return derivnum;
					case 1:
						return num * (a1 * R / ((a1 + a2) * (a1 + a2)) - R / (a1 + a2)) +
								derivnum * -R * a1 / (a1 + a2);
				}
			case 1:
				switch (l2) {
					case 0:
						return -num * R * a2 / ((a1 + a2) * (a1 + a2)) + derivnum * R * a2 / (a1 + a2);
					case 1:
						return (1 / (2 * (a1 + a2)) - a1 * a2 * R * R / ((a1 + a2) * (a1 + a2))) * derivnum
								+ num * (2 * a1 * a2 * R * R / ((a1 + a2) * (a1 + a2) * (a1 + a2)) -
								a2 * R * R / ((a1 + a2) * (a1 + a2)) - 1 / (2 * (a1 + a2) * (a1 + a2)));
				}
		}

		return 0;
	}

	private static double Ialphadiagp2d(int l1, int l2, double a1, double a2, double R) {
		double exp = Pow.exp(-a1 * a2 * R * R / (a1 + a2));
		double sqrt = Math.sqrt(Math.PI / (a1 + a2));
		double num = sqrt * exp;

		double derivnum = -exp * R * R * (a2 / (a1 + a2) - a1 * a2 / ((a1 + a2) * (a1 + a2))) * sqrt -
				exp * 0.5 * Math.sqrt(Math.PI) / Pow.pow(a1 + a2, 1.5);

		double derivnum2 = exp * Math.sqrt(Math.PI) * Pow.pow(a1 + a2, -4.5) *
				(Pow.pow(R, 4) * Pow.pow(a2, 4) + 3 * R * R * a2 * a2 * (a1 + a2) + 0.75 * (a1 + a2) * (a1 + a2));

		double v2 = (a1 + a2) * (a1 + a2) * (a1 + a2);
		double v = num * 2 * a2 * R / v2;
		double v1 = 2 * a2 * R / ((a1 + a2) * (a1 + a2)) * derivnum;
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return derivnum2;
					case 1:
						return v - v1 - a1 * R / (a1 + a2) * derivnum2;
				}
			case 1:
				switch (l2) {
					case 0:
						return v - v1 + a2 * R / (a1 + a2) * derivnum2;
					case 1:
						return (1 / (2 * (a1 + a2)) - a1 * a2 * R * R / ((a1 + a2) * (a1 + a2))) * derivnum2 +
								2 * derivnum * (2 * a1 * a2 * R * R / v2 - a2 * R * R / ((a1 + a2) * (a1 + a2)) -
										1 / (2 * (a1 + a2) * (a1 + a2))) +
								num * ((4 * a2 * a2 - 2 * a1 * a2) * R * R / Pow.pow(a1 + a2, 4) + Pow.pow(a1 + a2, -3));
				}
		}

		return 0;
	}

	private static double Ialphacrossp2d(int l1, int l2, double a1, double a2, double R) {
		double exp = Pow.exp(-a1 * a2 * R * R / (a1 + a2));
		double sqrt = Math.sqrt(Math.PI / (a1 + a2));
		double num = sqrt * exp;

		double v = a1 * a2 / ((a1 + a2) * (a1 + a2));
		double derivnuma = -exp * R * R * (a2 / (a1 + a2) - v) *
				sqrt - exp * 0.5 * Math.sqrt(Math.PI) / Pow.pow(a1 + a2, 1.5);

		double derivnumb = -exp * R * R * (a1 / (a1 + a2) - v) *
				sqrt - exp * 0.5 * Math.sqrt(Math.PI) / Pow.pow(a1 + a2, 1.5);

		double derivnum2 =
				exp * Math.sqrt(Math.PI) * Pow.pow(a1 + a2, -4.5) * 0.25 *
						(4 * a1 * a1 * a2 * a2 * R * R * R * R +
								2 * R * R * (a1 + a2) * (a1 * a1 + a2 * a2 - 4 * a1 * a2) + 3 * (a1 + a2) * (a1 + a2));

		double v3 = (a1 + a2) * (a1 + a2) * (a1 + a2);
		double v4 = R / ((a1 + a2) * (a1 + a2));
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return derivnum2;
					case 1:
						return (v4 - 2 * a1 * R / v3) * num
								+ (a1 * R / ((a1 + a2) * (a1 + a2)) - R / (a1 + a2)) * derivnumb
								+ R * a1 / ((a1 + a2) * (a1 + a2)) * derivnuma
								+ derivnum2 * -R * a1 / (a1 + a2);

				}
			case 1:
				switch (l2) {
					case 0:
						return -(v4 - 2 * a2 * R / v3) * num
								- (a2 * R / ((a1 + a2) * (a1 + a2)) - R / (a1 + a2)) * derivnuma
								- R * a2 / ((a1 + a2) * (a1 + a2)) * derivnumb
								- derivnum2 * -R * a2 / (a1 + a2);
					case 1:
						double v1 = 2 * a1 * a2 * R * R / v3;
						double v2 = 1 / (2 * (a1 + a2) * (a1 + a2));
						return (1 / (2 * (a1 + a2)) - a1 * a2 * R * R / ((a1 + a2) * (a1 + a2))) * derivnum2
								+ derivnumb * (v1 -
								a2 * R * R / ((a1 + a2) * (a1 + a2)) - v2)
								+ derivnuma * (v1 -
								a1 * R * R / ((a1 + a2) * (a1 + a2)) - v2)
								+ num * ((a1 * a1 - 4 * a1 * a2 + a2 * a2) * R * R /
								Pow.pow(a1 + a2, 4) + 1 / v3);
				}
		}

		return 0;
	}

	private static double Ialphapgd(int l1, int l2, double a1, double a2, double R) {
		double exp = Pow.exp(-a1 * a2 * R * R / (a1 + a2));
		double v = a2 / (a1 + a2) - a1 * a2 / ((a1 + a2) * (a1 + a2));
		double sqrt = Math.sqrt(Math.PI / (a1 + a2));
		double derivnumalpha = -exp * R * R * v * sqrt - exp * 0.5 * Math.sqrt(Math.PI) / Pow.pow(a1 + a2, 1.5);

		double num = sqrt * exp;

		double derivnum = num * 2 * a1 * a2 * R / (a1 + a2);

		double derivnum2 = derivnumalpha * R * 2 * a1 * a2 / (a1 + a2) +
				2 * v * R * num;

		double v2 = a2 / ((a1 + a2) * (a1 + a2));
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return derivnum2;
					case 1:
						return a1 / (a1 + a2) * (derivnumalpha - R * derivnum2) +
								v2 * (num - R * derivnum);
				}
			case 1:
				switch (l2) {
					case 0:
						return a2 / (a1 + a2) * (R * derivnum2 - derivnumalpha) -
								v2 * (R * derivnum - num);
					case 1:
						double v1 = (a1 + a2) * (a1 + a2) * (a1 + a2);
						return (1 / (2 * (a1 + a2)) - a1 * a2 * R * R / ((a1 + a2) * (a1 + a2))) * derivnum2
								+ derivnumalpha * 2 * a1 * a2 * R / ((a1 + a2) * (a1 + a2))
								+ num * (2 * a2 * R / ((a1 + a2) * (a1 + a2)) -
								4 * a1 * a2 * R / v1)
								+ derivnum * (2 * a1 * a2 * R * R / v1 -
								a2 * R * R / ((a1 + a2) * (a1 + a2)) - 1 / (2 * (a1 + a2) * (a1 + a2)));

				}
		}

		return 0;
	}

	private static double Ialphacrossp2gd(int l1, int l2, double a1, double a2, double R) {
		double exp = Pow.exp(-a1 * a2 * R * R / (a1 + a2));
		double v2 = a1 * a2 / ((a1 + a2) * (a1 + a2));
		double v = a2 / (a1 + a2) - v2;
		double sqrt = Math.sqrt(Math.PI / (a1 + a2));
		double derivnumalphaa = -exp * R * R * v * sqrt - exp * 0.5 * Math.sqrt(Math.PI) / Pow.pow(a1 + a2, 1.5);
		double v1 = a1 / (a1 + a2) - v2;
		double derivnumalphab = -exp * R * R * v1 * sqrt - exp * 0.5 * Math.sqrt(Math.PI) / Pow.pow(a1 + a2, 1.5);

		double num = sqrt * exp;

		double derivnumgeom = num * 2 * a1 * a2 * R / (a1 + a2);

		double derivnumalphageoma = derivnumalphaa * R * 2 * a1 * a2 / (a1 + a2) +
				2 * v * R * num;

		double derivnumalphageomb = derivnumalphab * R * 2 * a1 * a2 / (a1 + a2) +
				2 * v1 * R * num;

		double v4 = a1 * a1 + a2 * a2 - 4 * a1 * a2;
		double v3 = 4 * a1 * a1 * a2 * a2 * R * R * R * R +
				2 * R * R * (a1 + a2) * v4 + 3 * (a1 + a2) * (a1 + a2);
		double derivnumalpha2cross = exp * Math.sqrt(Math.PI) * Pow.pow(a1 + a2, -4.5) * 0.25 * v3;

		double derivnum3 = 0.25 * Pow.pow(a1 + a2, -4) * derivnumgeom * v3 -
				num * 0.25 * Pow.pow(a1 + a2, -4) * (16 * a1 * a1 * a2 * a2 * R * R * R + 4 * R * (a1 + a2) * v4);

		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return derivnum3;
					case 1:
						return (2 * a1 * Pow.pow(a1 + a2, -3) - Pow.pow(a1 + a2, -2)) * num +
								(R * Pow.pow(a1 + a2, -2) - 2 * a1 * R * Pow.pow(a1 + a2, -3)) * derivnumgeom
								+ (Pow.pow(a1 + a2, -1) - a1 * Pow.pow(a1 + a2, -2)) * derivnumalphab +
								(a1 * R * Pow.pow(a1 + a2, -2) - R * Pow.pow(a1 + a2, -1)) * derivnumalphageomb
								- a1 * Pow.pow(a1 + a2, -2) * derivnumalphaa +
								a1 * R * Pow.pow(a1 + a2, -2) * derivnumalphageoma
								+ a1 * Pow.pow(a1 + a2, -1) * derivnumalpha2cross -
								a1 * R * Pow.pow(a1 + a2, -1) * derivnum3;
				}
			case 1:
				switch (l2) {
					case 0:
						return (R / (a1 + a2) - a2 * R * Pow.pow(a1 + a2, -2)) * derivnumalphageoma -
								(1 / (a1 + a2) - a2 * Pow.pow(a1 + a2, -2)) * derivnumalphaa
								- a2 / (a1 + a2) * derivnumalpha2cross + a2 * R / (a1 + a2) * derivnum3
								+ (2 * a2 * R * Pow.pow(a1 + a2, -3) - R * Pow.pow(a1 + a2, -2)) * derivnumgeom -
								(2 * a2 * Pow.pow(a1 + a2, -3) - Pow.pow(a1 + a2, -2)) * num
								+ a2 * Pow.pow(a1 + a2, -2) * derivnumalphab -
								a2 * R * Pow.pow(a1 + a2, -2) * derivnumalphageomb;
					case 1:
						return (0.5 / (a1 + a2) - a1 * a2 * R * R * Pow.pow(a1 + a2, -2)) * derivnum3 +
								2 * a1 * a2 * R * Pow.pow(a1 + a2, -2) * derivnumalpha2cross
								+ derivnumalphaa *
								(2 * a1 * R * Pow.pow(a1 + a2, -2) - 4 * a1 * a2 * R * Pow.pow(a1 + a2, -3)) +
								derivnumalphageoma *
										(2 * a1 * a2 * R * R * Pow.pow(a1 + a2, -3) - 0.5 * Pow.pow(a1 + a2, -2) -
												a1 * R * R * Pow.pow(a1 + a2, -2))
								+ derivnumalphab *
								(2 * a2 * R * Pow.pow(a1 + a2, -2) - 4 * a1 * a2 * R * Pow.pow(a1 + a2, -3)) +
								derivnumalphageomb * (2 * a1 * a2 * R * R * Pow.pow(a1 + a2, -3) -
										a2 * R * R * Pow.pow(a1 + a2, -2) - 0.5 * Pow.pow(a1 + a2, -2))
								+ derivnumgeom * (v4 * R * R * Pow.pow(a1 + a2, -4) +
								Pow.pow(a1 + a2, -3))
								- 2 * v4 * R * Pow.pow(a1 + a2, -4) * num;

				}
		}

		return 0;
	}

	private static double Ialphadiagp2gd(int l1, int l2, double a1, double a2, double R) {

		double exp = Pow.exp(-a1 * a2 * R * R / (a1 + a2));
		double v = a2 / (a1 + a2) - a1 * a2 / ((a1 + a2) * (a1 + a2));
		double sqrt = Math.sqrt(Math.PI / (a1 + a2));
		double derivnumalpha = -exp * R * R * v * sqrt - exp * 0.5 * Math.sqrt(Math.PI) / Pow.pow(a1 + a2, 1.5);

		double num = sqrt * exp;

		double derivnumgeom = num * 2 * a1 * a2 * R / (a1 + a2);

		double derivnumalphageom = derivnumalpha * R * 2 * a1 * a2 / (a1 + a2) +
				2 * v * R * num;

		double v3 = 0.75 * (a1 + a2) * (a1 + a2);
		double v2 = 3 * R * R * a2 * a2 * (a1 + a2);
		double derivnumalpha2diag =
				exp * Math.sqrt(Math.PI) * Pow.pow(a1 + a2, -4.5) * (Pow.pow(R, 4) * Pow.pow(a2, 4) + v2 + v3);

		double derivnum3 = Pow.pow(a1 + a2, -4) * derivnumgeom * (a2 * a2 * a2 * a2 * R * R * R * R + v2 + v3) -
				num * Pow.pow(a1 + a2, -4) * (4 * a2 * a2 * a2 * a2 * R * R * R + 6 * R * a2 * a2 * (a1 + a2));

		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return derivnum3;
					case 1:
						return 2 * a2 * R * Pow.pow(a1 + a2, -3) * derivnumgeom -
								2 * a2 * Pow.pow(a1 + a2, -3) * num
								+ 2 * a2 * Pow.pow(a1 + a2, -2) * derivnumalpha -
								2 * a2 * R * Pow.pow(a1 + a2, -2) * derivnumalphageom
								+ a1 / (a1 + a2) * derivnumalpha2diag - a1 * R / (a1 + a2) * derivnum3;
				}
			case 1:
				switch (l2) {
					case 0:
						return 2 * a2 * R * Pow.pow(a1 + a2, -3) * derivnumgeom -
								2 * a2 * Pow.pow(a1 + a2, -3) * num
								+ 2 * a2 * Pow.pow(a1 + a2, -2) * derivnumalpha -
								2 * a2 * R * Pow.pow(a1 + a2, -2) * derivnumalphageom
								- a2 / (a1 + a2) * derivnumalpha2diag + a2 * R / (a1 + a2) * derivnum3;
					case 1:
						double v1 = 4 * a2 * a2 - 2 * a1 * a2;
						return (0.5 / (a1 + a2) - a1 * a2 * R * R * Pow.pow(a1 + a2, -2)) * derivnum3 +
								2 * a1 * a2 * R * Pow.pow(a1 + a2, -2) * derivnumalpha2diag
								+ 2 * derivnumalphageom *
								(2 * a1 * a2 * R * R * Pow.pow(a1 + a2, -3) - a2 * R * R * Pow.pow(a1 + a2, -2) -
										0.5 * Pow.pow(a1 + a2, -2))
								- 2 * derivnumalpha *
								(4 * a1 * a2 * R * Pow.pow(a1 + a2, -3) - 2 * a2 * R * Pow.pow(a1 + a2, -2))
								+ derivnumgeom *
								(v1 * R * R * Pow.pow(a1 + a2, -4) + Pow.pow(a1 + a2, -3))
								- 2 * v1 * R * Pow.pow(a1 + a2, -4) * num;

				}
		}

		return 0;
	}

	public double Nalphapd() {
		return (2 * L + 3) / 4.0 * Pow.pow(2 / Math.PI, 0.75) * Pow.pow(2, L) *
				Pow.pow(exponent, (2 * L - 1) / 4.0) /
				Math.sqrt(fact2(2 * i - 1) * fact2(2 * j - 1) * fact2(2 * k - 1));
	}

	public double Nalphap2d() {
		return (2 * L + 3) * (2 * L - 1) / 16.0 * Pow.pow(2 / Math.PI, 0.75) * Pow.pow(2, L) *
				Pow.pow(exponent, (2 * L - 5) / 4.0) /
				Math.sqrt(fact2(2 * i - 1) * fact2(2 * j - 1) * fact2(2 * k - 1));
	}

	public double getN() {
		return this.N;
	} // todo remove getter

	public boolean equals(Object b) {
		if (b instanceof GTO) {
			GTO c = (GTO) b;
			return this.i == c.i && this.j == c.j && this.k == c.k &&
					this.exponent == c.exponent &&
					Arrays.equals(this.coordinates, c.coordinates);
		}
		return false;
	}
}
