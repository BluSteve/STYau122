package nddo.scf;

import nddo.structs.OrbitalProperties;
import tools.Pow;

import java.util.Arrays;

// this is a complete orbital representation
public class LCGTO extends Orbital { // HAS AtomFixed, OrbitalProperties, c, e
	//LCGTO = Linear Combinations (of) GTOs i.e. contracted basis functions
	//Each LCGTO is described by an array of GTOs and an array of contraction
	// coefficients

	protected final double N; //normalisation coefficient, should be very close to 1.
	protected final int n; //number of Gaussians in this contraction
	protected final GTO[] gaussArray;
	protected final double[] gaussExponents;
	protected final double[] coefficientArray;

	public LCGTO(OrbitalProperties op, double[] coordinates, double[] e, double[] c) {
		// e = exponent array, c = coefficient array
		super(op, coordinates);
		this.gaussExponents = e;
		this.coefficientArray = c;
		this.n = c.length;

		gaussArray = new GTO[n];

		for (int a = 0; a < n; a++) {
			GTO g = new GTO(this.op, this.coordinates, e[a]);
			gaussArray[a] = g;
		}

		double sum = 0;

		for (int a = 0; a < n; a++) {
			for (int b = 0; b < n; b++) {
				sum += c[a] * c[b] * gaussArray[a].getN() * gaussArray[b].getN() / Pow.pow(e[a] + e[b], L + 1.5);
			}
		}

		sum *= Pow.pow(Math.PI, 1.5) * GTO.fact2(2 * i - 1) * GTO.fact2(2 * j - 1) * GTO.fact2(2 * k - 1)
				/ Pow.pow(2, L);
		sum = Pow.pow(sum, -0.5);
		N = sum;
	}

	public LCGTO(LCGTO lcgto) {
		super(lcgto);

		this.gaussExponents = lcgto.gaussExponents;
		this.coefficientArray = lcgto.coefficientArray;
		this.n = lcgto.n;
		this.gaussArray = lcgto.gaussArray;
		this.N = lcgto.N;
	}

	public static double S(LCGTO X1, LCGTO X2) {//normalised overlap integral
		double S = 0;

		for (int i = 0; i < X1.getn(); i++) {
			for (int j = 0; j < X2.getn(); j++) {
				S += X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
						GTO.S(X1.getGaussArray()[i], X2.getGaussArray()[j]);
			}

		}

		return S * X1.getN() * X2.getN();
	}

	public static double Sgd(LCGTO X1, LCGTO X2, int tau) {
		double S = 0;

		for (int i = 0; i < X1.getn(); i++) {
			for (int j = 0; j < X2.getn(); j++) {
				S += X1.getCoeffArray()[i] * X2.getCoeffArray()[j] * GTO.Sgd(X1.getGaussArray()[i], X2.getGaussArray()[j], tau);
			}
		}

		return S * X1.getN() * X2.getN();
	}

	public static double Sg2d(LCGTO X1, LCGTO X2, int tau1, int tau2) {
		double S = 0;

		for (int i = 0; i < X1.getn(); i++) {
			for (int j = 0; j < X2.getn(); j++) {
				S += X1.getCoeffArray()[i] * X2.getCoeffArray()[j] * GTO.Sg2d(X1.getGaussArray()[i], X2.getGaussArray()[j], tau1, tau2);
			}

		}

		return S * X1.getN() * X2.getN();
	}

	public GTO[] getGaussArray() {
		return this.gaussArray;
	}

	public int getn() {
		return n;
	}

	public double[] getCoeffArray() {
		return coefficientArray;
	}

	public double[] getCoords() {
		return this.coordinates;
	}

	public int getL() {
		return L;
	}

	public int geti() {
		return i;
	}

	public int getj() {
		return this.j;
	}

	public int getk() {
		return this.k;
	}

	public double getN() {
		return this.N;
	}

	public double[] getGaussExponents() {
		return gaussExponents;
	}

	public boolean equals(Object b) {
		if (b instanceof LCGTO) {
			LCGTO a = (LCGTO) b;
			return this.N == a.N && Arrays.equals(a.coordinates,
					coordinates) &&
					a.i == i && a.j == j && a.k == k &&
					Arrays.equals(a.coefficientArray, coefficientArray) &&
					Arrays.deepEquals(a.gaussArray, gaussArray);

		}
		return false;
	}
}
