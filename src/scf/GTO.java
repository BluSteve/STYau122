package scf;

import java.util.Arrays;

//this is the GTO class, and all the relevant GTO integral routines are
// implemented here.
//The GTO object represents a GTO of the form (x-coordinates[0])^i
// (y-coordinates[1])^j
// (z-coordinates[2])^k exp(-exponent * r^2)
public class GTO {

	private int i, j, k;// angular terms
	private double N;
	private int L; // L = i + j + k
	private double exponent; //radial exponent
	private double[] coordinates; //coordinates

	public GTO(int i, int j, int k, double exponent, double[] coordinates) {
		this.i = i;
		this.j = j;
		this.k = k;
		this.L = i + j + k;
		this.exponent = exponent;
		this.coordinates = coordinates;

		this.N = Math.pow(2 / Math.PI, 0.75) * Math.pow(2, L) *
				Math.pow(exponent, (2 * L + 3) / 4.0) /
				Math.sqrt(
						fact2(2 * i - 1) * fact2(2 * j - 1) * fact2(2 * k - 1));
	}

	static int fact2(int num) {
		int sum = 1;
		for (int i = num; i > 0; i -= 2) {
			sum *= i;
		}

		return sum;
	}

	public static double R(double[] P, double[] C) {

		double val =
				(P[0] - C[0]) * (P[0] - C[0]) + (P[1] - C[1]) * (P[1] - C[1]) +
						(P[2] - C[2]) * (P[2] - C[2]);

		return Math.sqrt(val);


	}

	public static double getS(GTO a, GTO b) {

		return a.getN() * b.getN() *
				I(a.i, b.i, a.exponent, b.exponent,
						b.coordinates[0] - a.coordinates[0]) *
				I(a.j, b.j, a.exponent, b.exponent,
						b.coordinates[1] - a.coordinates[1]) *
				I(a.k, b.k, a.exponent, b.exponent,
						b.coordinates[2] - a.coordinates[2]);
	}

	public static double getSderiv(GTO a, GTO b, int tau) {

		switch (tau) {
			case 0:

				return a.getN() * b.getN() *
						Ideriv(a.i, b.i, a.exponent, b.exponent,
								b.coordinates[0] - a.coordinates[0]) *
						I(a.j, b.j, a.exponent, b.exponent,
								b.coordinates[1] - a.coordinates[1]) *
						I(a.k, b.k, a.exponent, b.exponent,
								b.coordinates[2] - a.coordinates[2]);
			case 1:
				return a.getN() * b.getN() * I(a.i, b.i, a.exponent,
						b.exponent,
						b.coordinates[0] - a.coordinates[0]) *
						Ideriv(a.j, b.j, a.exponent, b.exponent,
								b.coordinates[1] - a.coordinates[1]) *
						I(a.k, b.k, a.exponent, b.exponent,
								b.coordinates[2] - a.coordinates[2]);
			case 2:
				return a.getN() * b.getN() * I(a.i, b.i, a.exponent,
						b.exponent,
						b.coordinates[0] - a.coordinates[0]) *
						I(a.j, b.j, a.exponent, b.exponent,
								b.coordinates[1] - a.coordinates[1]) *
						Ideriv(a.k, b.k, a.exponent, b.exponent,
								b.coordinates[2] - a.coordinates[2]);
		}
		return 0;
	}

	public static double getSderiv2(GTO a, GTO b, int tau1, int tau2) {

		int A = Math.min(tau1, tau2);
		int B = Math.max(tau1, tau2);

		switch (A) {
			case 0:
				switch (B) {
					case 0: //derivative wrt x and x
						return a.getN() * b.getN() *
								Ideriv2(a.i, b.i, a.exponent, b.exponent,
										b.coordinates[0] - a.coordinates[0]) *
								I(a.j, b.j, a.exponent, b.exponent,
										b.coordinates[1] - a.coordinates[1]) *
								I(a.k, b.k, a.exponent, b.exponent,
										b.coordinates[2] - a.coordinates[2]);
					case 1: // derivative wrt x and y
						return a.getN() * b.getN() *
								Ideriv(a.i, b.i, a.exponent, b.exponent,
										b.coordinates[0] - a.coordinates[0]) *
								Ideriv(a.j, b.j, a.exponent, b.exponent,
										b.coordinates[1] - a.coordinates[1]) *
								I(a.k, b.k, a.exponent, b.exponent,
										b.coordinates[2] - a.coordinates[2]);
					case 2: // derivative wrt x and z
						return a.getN() * b.getN() *
								Ideriv(a.i, b.i, a.exponent, b.exponent,
										b.coordinates[0] - a.coordinates[0]) *
								I(a.j, b.j, a.exponent, b.exponent,
										b.coordinates[1] - a.coordinates[1]) *
								Ideriv(a.k, b.k, a.exponent, b.exponent,
										b.coordinates[2] - a.coordinates[2]);

				}
			case 1:
				switch (B) {
					case 1: // derivative wrt y and y
						return a.getN() * b.getN() *
								I(a.i, b.i, a.exponent, b.exponent,
										b.coordinates[0] - a.coordinates[0]) *
								Ideriv2(a.j, b.j, a.exponent, b.exponent,
										b.coordinates[1] - a.coordinates[1]) *
								I(a.k, b.k, a.exponent, b.exponent,
										b.coordinates[2] - a.coordinates[2]);
					case 2: // derivative wrt y and z
						return a.getN() * b.getN() *
								I(a.i, b.i, a.exponent, b.exponent,
										b.coordinates[0] - a.coordinates[0]) *
								Ideriv(a.j, b.j, a.exponent, b.exponent,
										b.coordinates[1] - a.coordinates[1]) *
								Ideriv(a.k, b.k, a.exponent, b.exponent,
										b.coordinates[2] - a.coordinates[2]);

				}
			case 2: // derivative wrt z and z
				return a.getN() * b.getN() * I(a.i, b.i, a.exponent,
						b.exponent,
						b.coordinates[0] - a.coordinates[0]) *
						I(a.j, b.j, a.exponent, b.exponent,
								b.coordinates[1] - a.coordinates[1]) *
						Ideriv2(a.k, b.k, a.exponent, b.exponent,
								b.coordinates[2] - a.coordinates[2]);
		}

		return 0;

	}

	private static double I(int l1, int l2, double a1, double a2, double R) {

		double num =
				Math.sqrt(Math.PI / (a1 + a2)) *
						Math.exp(-a1 * a2 * R * R / (a1 + a2));
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
						return (1 / (2 * (a1 + a2)) -
								a1 * a2 * R * R / ((a1 + a2) * (a1 + a2))) *
								num;
				}
				break;
		}

		return 0;
	}

	private static double Ideriv(int l1, int l2, double a1, double a2,
								 double R) {

		double num =
				Math.sqrt(Math.PI / (a1 + a2)) *
						Math.exp(-a1 * a2 * R * R / (a1 + a2));
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return num * 2 * a1 * a2 * R / (a1 + a2);
					case 1:
						return a1 / (a1 + a2) * num *
								(1 - 2 * a1 * a2 * R * R / (a1 + a2));
				}
				break;
			case 1:
				switch (l2) {
					case 0:
						return a2 / (a1 + a2) * num *
								(2 * a1 * a2 * R * R / (a1 + a2) - 1);
					case 1:
						return num * a1 * a2 * R / Math.pow(a1 + a2, 2) *
								(3 - 2 * R * R * a1 * a2 / (a1 + a2));
				}
				break;
		}

		return 0;

	}

	private static double Ideriv2(int l1, int l2, double a1, double a2,
								  double R) {

		double num =
				Math.sqrt(Math.PI / (a1 + a2)) *
						Math.exp(-a1 * a2 * R * R / (a1 + a2));
		switch (l1) {
			case 0:
				switch (l2) {
					case 0:
						return num * 2 * a1 * a2 / (a1 + a2) *
								(2 * a1 * a2 / (a1 + a2) * R * R - 1);
					case 1:
						return 2 * a1 * a1 * a2 * R / ((a1 + a2) * (a1 + a2)) *
								num *
								(3 - 2 * a1 * a2 / (a1 + a2) * R * R);
				}
				break;
			case 1:
				switch (l2) {
					case 0:
						return 2 * a1 * a2 * a2 * R / ((a1 + a2) * (a1 + a2)) *
								num *
								(2 * a1 * a2 / (a1 + a2) * R * R - 3);
					case 1:
						return a1 * a2 * Math.pow(a1 + a2, -2) * num *
								(2 * R * a1 * a2 / (a1 + a2) *
										(3 * R - 2 * R * R * R * a1 * a2 /
												(a1 + a2)) +
										6 * R * R * a1 * a2 / (a1 + a2) - 3);
				}
				break;
		}

		return 0;

	}

	public static double getSderivfinite(GTO a, GTO b, int tau) {

		double[] newcoords = a.coordinates.clone();
		newcoords[tau] += 1E-8;

		GTO newa = new GTO(a.i, a.j, a.k, a.exponent, newcoords);

		return 1E8 * (GTO.getS(newa, b) - GTO.getS(a, b));

	}

	public static double getSderiv2finite(GTO a, GTO b, int tau1, int tau2) {
		double[] newcoords = a.coordinates.clone();
		newcoords[tau2] += 1E-8;

		GTO newa = new GTO(a.i, a.j, a.k, a.exponent, newcoords);

		return 1E8 * (GTO.getSderiv(newa, b, tau1) - GTO.getSderiv(a, b,
				tau1));
	}

	public int geti() {
		return i;
	}

	public int getj() {
		return j;
	}

	public int getk() {
		return k;
	}

	public double getexponent() {
		return exponent;
	}

	public double getX() {
		return coordinates[0];
	}

	public double getY() {
		return coordinates[1];
	}

	public double getZ() {
		return coordinates[2];
	}

	public double[] getcoords() {
		return this.coordinates;
	}

	public double getN() {
		return this.N;
	}

	public int getL() {
		return this.L;
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
