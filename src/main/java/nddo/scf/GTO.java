package nddo.scf;

import nddo.structs.OrbitalProperties;

import java.util.Arrays;

/*
this is the GTO class, and all the relevant GTO integral routines are implemented here. The GTO object represents a
GTO of the form (x-coordinates[0])^i (y-coordinates[1])^j (z-coordinates[2])^k exp(-exponent * r^2)*/
public class GTO extends Orbital {
	private final double N;
	private final double exponent; /*radial exponent*/

	public GTO(OrbitalProperties op, double[] coordinates, double exponent) {
		super(op, coordinates);
		this.exponent = exponent;

		this.N = Math.pow(2 / Math.PI, 0.75) * Math.pow(2, L) *
				Math.pow(exponent, (2 * L + 3) / 4.0) /
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

	public static double getS(GTO a, GTO b) {
		return a.getN() * b.getN() * I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
				I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
				I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
	}

	public static double getSderiv(GTO a, GTO b, int tau) {
		switch (tau) {
			case 0:
				return a.getN() * b.getN() *
						Ideriv(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
						I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
						I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
			case 1:
				return a.getN() * b.getN() * I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
						Ideriv(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
						I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
			case 2:
				return a.getN() * b.getN() * I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
						I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
						Ideriv(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
		}
		return 0;
	}

	public static double getSderiv2(GTO a, GTO b, int tau1, int tau2) {
		int A = Math.min(tau1, tau2);
		int B = Math.max(tau1, tau2);

		switch (A) {
			case 0:
				switch (B) {
					case 0: /*derivative wrt x and x*/
						return a.getN() * b.getN() *
								Ideriv2(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
					case 1: /* derivative wrt x and y*/
						return a.getN() * b.getN() *
								Ideriv(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								Ideriv(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
					case 2: /* derivative wrt x and z*/
						return a.getN() * b.getN() *
								Ideriv(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								Ideriv(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

				}
			case 1:
				switch (B) {
					case 1: /* derivative wrt y and y*/
						return a.getN() * b.getN() *
								I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								Ideriv2(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
					case 2: /* derivative wrt y and z*/
						return a.getN() * b.getN() *
								I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
								Ideriv(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
								Ideriv(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
				}
			case 2: /* derivative wrt z and z*/
				return a.getN() * b.getN() * I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0]) *
						I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1]) *
						Ideriv2(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
		}
		return 0;
	}

	public static double getSderivalpha(GTO a, GTO b, int type) {
		switch (type) {
			case 0:
				return a.Nderivalpha() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Iderivalpha(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Iderivalpha(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Iderivalpha(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

			case 1:
				return a.getN() * b.Nderivalpha()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Iderivalpha(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Iderivalpha(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Iderivalpha(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2]);
		}

		return 0;
	}

	public static double getSderiv2alphadiag(GTO a, GTO b, int type) {
		switch (type) {
			case 0:
				return a.Nderiv2alpha() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.Nderivalpha() * b.getN()
						* Iderivalpha(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.Nderivalpha() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Iderivalpha(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.Nderivalpha() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Iderivalpha(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* Iderivalpha(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Iderivalpha(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* Iderivalpha(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Iderivalpha(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Iderivalpha(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Iderivalpha(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* Ideriv2alphadiag(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ideriv2alphadiag(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ideriv2alphadiag(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);

			case 1:
				return a.getN() * b.Nderiv2alpha()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.Nderivalpha()
						* Iderivalpha(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.Nderivalpha()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Iderivalpha(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.Nderivalpha()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Iderivalpha(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* Iderivalpha(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* Iderivalpha(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* Iderivalpha(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Iderivalpha(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ 2 * a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Iderivalpha(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* Iderivalpha(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

						+ a.getN() * b.getN()
						* Ideriv2alphadiag(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* Ideriv2alphadiag(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
						* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

						+ a.getN() * b.getN()
						* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
						* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
						* Ideriv2alphadiag(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2]);
		}

		return 0;
	}

	public static double getSderiv2alphacross(GTO a, GTO b) {
		return a.Nderivalpha() * b.Nderivalpha()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.Nderivalpha()
				* Iderivalpha(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.Nderivalpha()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Iderivalpha(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.Nderivalpha()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Iderivalpha(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.Nderivalpha() * b.getN()
				* Iderivalpha(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.Nderivalpha() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Iderivalpha(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.Nderivalpha() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Iderivalpha(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

				+ a.getN() * b.getN()
				* Iderivalpha(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Iderivalpha(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* Iderivalpha(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Iderivalpha(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

				+ a.getN() * b.getN()
				* Iderivalpha(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
				* Iderivalpha(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Iderivalpha(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Iderivalpha(b.k, a.k, b.exponent, a.exponent, a.coordinates[2] - b.coordinates[2])

				+ a.getN() * b.getN()
				* Iderivalpha(b.i, a.i, b.exponent, a.exponent, a.coordinates[0] - b.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Iderivalpha(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Iderivalpha(b.j, a.j, b.exponent, a.exponent, a.coordinates[1] - b.coordinates[1])
				* Iderivalpha(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* Ideriv2alphacross(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* Ideriv2alphacross(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* I(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2])

				+ a.getN() * b.getN()
				* I(a.i, b.i, a.exponent, b.exponent, b.coordinates[0] - a.coordinates[0])
				* I(a.j, b.j, a.exponent, b.exponent, b.coordinates[1] - a.coordinates[1])
				* Ideriv2alphacross(a.k, b.k, a.exponent, b.exponent, b.coordinates[2] - a.coordinates[2]);
	}

	private static double I(int l1, int l2, double a1, double a2, double R) {
		double num = Math.sqrt(Math.PI / (a1 + a2)) * Math.exp(-a1 * a2 * R * R / (a1 + a2));
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

	private static double Ideriv(int l1, int l2, double a1, double a2, double R) {
		double num = Math.sqrt(Math.PI / (a1 + a2)) * Math.exp(-a1 * a2 * R * R / (a1 + a2));
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return num * 2 * a1 * a2 * R / (a1 + a2);
					case 1:
						return a1 / (a1 + a2) * num * (1 - 2 * a1 * a2 * R * R / (a1 + a2));
				}
				break;
			case 1:
				switch (l2) {
					case 0:
						return a2 / (a1 + a2) * num * (2 * a1 * a2 * R * R / (a1 + a2) - 1);
					case 1:
						return num * a1 * a2 * R / Math.pow(a1 + a2, 2) * (3 - 2 * R * R * a1 * a2 / (a1 + a2));
				}
				break;
		}
		return 0;
	}

	private static double Ideriv2(int l1, int l2, double a1, double a2, double R) {
		double num = Math.sqrt(Math.PI / (a1 + a2)) * Math.exp(-a1 * a2 * R * R / (a1 + a2));
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return num * 2 * a1 * a2 / (a1 + a2) * (2 * a1 * a2 / (a1 + a2) * R * R - 1);
					case 1:
						return 2 * a1 * a1 * a2 * R / ((a1 + a2) * (a1 + a2)) * num *
								(3 - 2 * a1 * a2 / (a1 + a2) * R * R);
				}
				break;
			case 1:
				switch (l2) {
					case 0:
						return 2 * a1 * a2 * a2 * R / ((a1 + a2) * (a1 + a2)) * num *
								(2 * a1 * a2 / (a1 + a2) * R * R - 3);
					case 1:
						return a1 * a2 * Math.pow(a1 + a2, -2) * num *
								(2 * R * a1 * a2 / (a1 + a2) * (3 * R - 2 * R * R * R * a1 * a2 / (a1 + a2)) +
										6 * R * R * a1 * a2 / (a1 + a2) - 3);
				}
				break;
		}
		return 0;
	}

	private static double Iderivalpha(int l1, int l2, double a1, double a2, double R) {
		double exp = Math.exp(-a1 * a2 * R * R / (a1 + a2));
		double sqrt = Math.sqrt(Math.PI / (a1 + a2));
		double derivnum = -exp * R * R *
				(a2 / (a1 + a2) - a1 * a2 / ((a1 + a2) * (a1 + a2))) * sqrt
				- exp * 0.5 * Math.sqrt(Math.PI) / Math.pow(a1 + a2, 1.5);

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

	private static double Ideriv2alphadiag(int l1, int l2, double a1, double a2, double R) {
		double exp = Math.exp(-a1 * a2 * R * R / (a1 + a2));
		double sqrt = Math.sqrt(Math.PI / (a1 + a2));
		double num = sqrt * exp;

		double derivnum = -exp * R * R *
				(a2 / (a1 + a2) - a1 * a2 / ((a1 + a2) * (a1 + a2))) * sqrt
				- exp * 0.5 * Math.sqrt(Math.PI) / Math.pow(a1 + a2, 1.5);

		double derivnum2 = exp * Math.sqrt(Math.PI) * Math.pow(a1 + a2, -4.5) *
				(Math.pow(R, 4) * Math.pow(a2, 4) + 3 * R * R * a2 * a2 * (a1 + a2) + 0.75 * (a1 + a2) * (a1 + a2));
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return derivnum2;
					case 1:
						return num * 2 * a2 * R / ((a1 + a2) * (a1 + a2) * (a1 + a2))
								- 2 * a2 * R / ((a1 + a2) * (a1 + a2)) * derivnum
								- a1 * R / (a1 + a2) * derivnum2;
				}
			case 1:
				switch (l2) {
					case 0:
						return num * 2 * a2 * R / ((a1 + a2) * (a1 + a2) * (a1 + a2))
								- 2 * a2 * R / ((a1 + a2) * (a1 + a2)) * derivnum
								+ a2 * R / (a1 + a2) * derivnum2;

					case 1:
						return (1 / (2 * (a1 + a2)) - a1 * a2 * R * R / ((a1 + a2) * (a1 + a2))) * derivnum2
								+ 2 * derivnum * (2 * a1 * a2 * R * R / ((a1 + a2) * (a1 + a2) * (a1 + a2)) -
								a2 * R * R / ((a1 + a2) * (a1 + a2)) - 1 / (2 * (a1 + a2) * (a1 + a2)))
								+ num *
								((4 * a2 * a2 - 2 * a1 * a2) * R * R / Math.pow(a1 + a2, 4) + Math.pow(a1 + a2, -3));


				}
		}

		return 0;
	}

	private static double Ideriv2alphacross(int l1, int l2, double a1, double a2, double R) {
		double exp = Math.exp(-a1 * a2 * R * R / (a1 + a2));
		double sqrt = Math.sqrt(Math.PI / (a1 + a2));
		double num = sqrt * exp;

		double derivnuma = -exp * R * R * (a2 / (a1 + a2) - a1 * a2 / ((a1 + a2) * (a1 + a2))) *
						sqrt - exp * 0.5 * Math.sqrt(Math.PI) / Math.pow(a1 + a2, 1.5);

		double derivnumb = -exp * R * R * (a1 / (a1 + a2) - a1 * a2 / ((a1 + a2) * (a1 + a2))) *
						sqrt - exp * 0.5 * Math.sqrt(Math.PI) / Math.pow(a1 + a2, 1.5);

		double derivnum2 =
				exp * Math.sqrt(Math.PI) * Math.pow(a1 + a2, -4.5) * 0.25 *
						(4 * a1 * a1 * a2 * a2 * R * R * R * R +
								2 * R * R * (a1 + a2) * (a1 * a1 + a2 * a2 - 4 * a1 * a2) + 3 * (a1 + a2) * (a1 + a2));

		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return derivnum2;
					case 1:
						return (R / ((a1 + a2) * (a1 + a2)) - 2 * a1 * R / ((a1 + a2) * (a1 + a2) * (a1 + a2))) * num
								+ (a1 * R / ((a1 + a2) * (a1 + a2)) - R / (a1 + a2)) * derivnumb
								+ R * a1 / ((a1 + a2) * (a1 + a2)) * derivnuma
								+ derivnum2 * -R * a1 / (a1 + a2);

				}
			case 1:
				switch (l2) {
					case 0:
						return -(R / ((a1 + a2) * (a1 + a2)) - 2 * a2 * R / ((a1 + a2) * (a1 + a2) * (a1 + a2))) * num
								- (a2 * R / ((a1 + a2) * (a1 + a2)) - R / (a1 + a2)) * derivnuma
								- R * a2 / ((a1 + a2) * (a1 + a2)) * derivnumb
								- derivnum2 * -R * a2 / (a1 + a2);
					case 1:
						return (1 / (2 * (a1 + a2)) - a1 * a2 * R * R / ((a1 + a2) * (a1 + a2))) * derivnum2
								+ derivnumb * (2 * a1 * a2 * R * R / ((a1 + a2) * (a1 + a2) * (a1 + a2)) -
								a2 * R * R / ((a1 + a2) * (a1 + a2)) - 1 / (2 * (a1 + a2) * (a1 + a2)))
								+ derivnuma * (2 * a1 * a2 * R * R / ((a1 + a2) * (a1 + a2) * (a1 + a2)) -
								a1 * R * R / ((a1 + a2) * (a1 + a2)) - 1 / (2 * (a1 + a2) * (a1 + a2)))
								+ num * ((a1 * a1 - 4 * a1 * a2 + a2 * a2) * R * R /
								((a1 + a2) * (a1 + a2) * (a1 + a2) * (a1 + a2)) +
								1 / ((a1 + a2) * (a1 + a2) * (a1 + a2)));

				}
		}

		return 0;
	}

	public double Nderivalpha() {
		return (2 * L + 3) / 4.0 * Math.pow(2 / Math.PI, 0.75) * Math.pow(2, L) *
				Math.pow(exponent, (2 * L - 1) / 4.0) /
				Math.sqrt(fact2(2 * i - 1) * fact2(2 * j - 1) * fact2(2 * k - 1));
	}

	public double Nderiv2alpha() {
		return (2 * L + 3) * (2 * L - 1) / 16.0 * Math.pow(2 / Math.PI, 0.75) * Math.pow(2, L) *
				Math.pow(exponent, (2 * L - 5) / 4.0) /
				Math.sqrt(fact2(2 * i - 1) * fact2(2 * j - 1) * fact2(2 * k - 1));
	}

	public double getExponent() {
		return exponent;
	}

	public double getN() {
		return this.N;
	}

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
