package nddo.math;

import nddo.scf.GTO;

public class Multipoles {
	public static double qq(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							double p22, double D12, double D22, double R) {
		double a00 = p01 + p02;
		return Math.pow(R * R + a00 * a00, -0.5);
	}

	public static double quz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							 double p22, double D12, double D22, double R) {
		double a01 = p01 + p12;
		return 0.5 * Math.pow((R + D12) * (R + D12) + a01 * a01, -0.5) -
				0.5 * Math.pow((R - D12) * (R - D12) + a01 * a01, -0.5);
	}

	public static double qQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		double a02 = p01 + p22;
		return 0.5 * Math.pow(R * R + 4 * D22 * D22 + a02 * a02, -0.5) - 0.5 * Math.pow(R * R + a02 * a02, -0.5);
	}

	public static double qQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							  double p22, double D12, double D22, double R) {
		double a02 = p01 + p22;
		return 0.25 * Math.pow((R + 2 * D22) * (R + 2 * D22) + a02 * a02, -0.5) +
				0.25 * Math.pow((R - 2 * D22) * (R - 2 * D22) + a02 * a02, -0.5) -
				0.5 * Math.pow(R * R + a02 * a02, -0.5);
	}

	public static double upiupi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		double a11 = p11 + p12;
		return 0.5 * Math.pow(R * R + (D11 - D12) * (D11 - D12) + a11 * a11, -0.5) -
				0.5 * Math.pow(R * R + (D11 + D12) * (D11 + D12) + a11 * a11, -0.5);
	}

	public static double uzuz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							  double p22, double D12, double D22, double R) {
		double a11 = p11 + p12;
		return 0.25 * Math.pow((R + D11 - D12) * (R + D11 - D12) + a11 * a11, -0.5) -
				0.25 * Math.pow((R + D11 + D12) * (R + D11 + D12) + a11 * a11, -0.5) -
				0.25 * Math.pow((R - D11 - D12) * (R - D11 - D12) + a11 * a11, -0.5) +
				0.25 * Math.pow((R - D11 + D12) * (R - D11 + D12) + a11 * a11, -0.5);
	}

	public static double upiQpiz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R) {
		double a12 = p11 + p22;
		return -0.25 * Math.pow((R - D22) * (R - D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -0.5) +
				0.25 * Math.pow((R - D22) * (R - D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -0.5) +
				0.25 * Math.pow((R + D22) * (R + D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -0.5) -
				0.25 * Math.pow((R + D22) * (R + D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -0.5);
	}

	public static double uzQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R) {
		double a12 = p11 + p22;
		return -0.25 * Math.pow((R + D11) * (R + D11) + 4 * D22 * D22 + a12 * a12, -0.5) +
				0.25 * Math.pow((R - D11) * (R - D11) + 4 * D22 * D22 + a12 * a12, -0.5) +
				0.25 * Math.pow((R + D11) * (R + D11) + a12 * a12, -0.5) -
				0.25 * Math.pow((R - D11) * (R - D11) + a12 * a12, -0.5);
	}

	public static double uzQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							   double p22, double D12, double D22, double R) {
		double a12 = p11 + p22;
		return -0.125 * Math.pow((R + D11 - 2 * D22) * (R + D11 - 2 * D22) + a12 * a12, -0.5) +
				0.125 * Math.pow((R - D11 - 2 * D22) * (R - D11 - 2 * D22) + a12 * a12, -0.5) -
				0.125 * Math.pow((R + D11 + 2 * D22) * (R + D11 + 2 * D22) + a12 * a12, -0.5) +
				0.125 * Math.pow((R - D11 + 2 * D22) * (R - D11 + 2 * D22) + a12 * a12, -0.5) +
				0.25 * Math.pow((R + D11) * (R + D11) + a12 * a12, -0.5) -
				0.25 * Math.pow((R - D11) * (R - D11) + a12 * a12, -0.5);
	}

	public static double QpipiQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R) {
		double a22 = p21 + p22;
		return 0.125 * Math.pow(R * R + 4 * (D21 - D22) * (D21 - D22) + a22 * a22, -0.5) +
				0.125 * Math.pow(R * R + 4 * (D21 + D22) * (D21 + D22) + a22 * a22, -0.5) -
				0.25 * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -0.5) -
				0.25 * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -0.5) + 0.25 * Math.pow(R * R + a22 * a22, -0.5);
	}

	public static double QxxQyy(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		double a22 = p21 + p22;
		return 0.25 * Math.pow(R * R + 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, -0.5) -
				0.25 * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -0.5) -
				0.25 * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -0.5) + 0.25 * Math.pow(R * R + a22 * a22, -0.5);
	}

	public static double QpipiQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		double a22 = p21 + p22;
		return 0.125 * Math.pow((R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 + a22 * a22, -0.5) +
				0.125 * Math.pow((R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 + a22 * a22, -0.5) -
				0.125 * Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -0.5) -
				0.125 * Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -0.5) -
				0.25 * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -0.5) + 0.25 * Math.pow(R * R + a22 * a22, -0.5);
	}

	public static double QzzQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		double a22 = p21 + p22;
		return 0.0625 * Math.pow((R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) + a22 * a22, -0.5) +
				0.0625 * Math.pow((R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) + a22 * a22, -0.5) +
				0.0625 * Math.pow((R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) + a22 * a22, -0.5) +
				0.0625 * Math.pow((R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) + a22 * a22, -0.5) -
				0.125 * Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -0.5) -
				0.125 * Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -0.5) -
				0.125 * Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -0.5) -
				0.125 * Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -0.5) +
				0.25 * Math.pow(R * R + a22 * a22, -0.5);
	}

	public static double QpizQpiz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		double a22 = p21 + p22;
		return 0.125 * Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -0.5) -
				0.125 * Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -0.5) -
				0.125 * Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -0.5) +
				0.125 * Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -0.5) -
				0.125 * Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -0.5) +
				0.125 * Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -0.5) +
				0.125 * Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -0.5) -
				0.125 * Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -0.5);
	}

	public static double ssss(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							  double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double ssppippi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double sspzpz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double ppippiss(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	public static double pzpzss(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	public static double ppippippippi(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QpipiQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double pxpxpypy(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QxxQyy(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double ppippipzpz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QpipiQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double pzpzppippi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QpipiQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	public static double pzpzpzpz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QzzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double spzss(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							   double p22, double D12, double D22, double R) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	public static double spzppippi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				uzQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double spzpzpz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				uzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double ssspz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							   double p22, double D12, double D22, double R) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double ppippispz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) -
				uzQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	public static double pzpzspz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) -
				uzQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	public static double sppisppi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		return upiupi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double spzspz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		return uzuz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double sppippipz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R) {
		return upiQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double ppipzsppi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R) {
		return -upiQpiz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	public static double ppipzppipz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R) {
		return QpizQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	public static double pxpypxpy(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		return 0.5 * (ppippippippi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) -
				pxpxpypy(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R));
	}

	public static double qqderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		double R = GTO.R(xA, xB);
		double a00 = p01 + p02;
		return (xB[tau] - xA[tau]) * Math.pow(R * R + a00 * a00, -1.5);
	}

	public static double quzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		double R = GTO.R(xA, xB);
		double a01 = p01 + p12;
		return 0.5 / R * (xB[tau] - xA[tau]) * (R + D12) * Math.pow((R + D12) * (R + D12) + a01 * a01, -1.5) -
				0.5 / R * (xB[tau] - xA[tau]) * (R - D12) * Math.pow((R - D12) * (R - D12) + a01 * a01, -1.5);
	}

	public static double qQpipideriv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		double R = GTO.R(xA, xB);
		double a02 = p01 + p22;
		return 0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D22 * D22 + a02 * a02, -1.5) -
				0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + a02 * a02, -1.5);
	}

	public static double qQzzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		double R = GTO.R(xA, xB);
		double a02 = p01 + p22;
		return 0.25 / R * (xB[tau] - xA[tau]) * (R + 2 * D22) *
				Math.pow((R + 2 * D22) * (R + 2 * D22) + a02 * a02, -1.5) +
				0.25 / R * (xB[tau] - xA[tau]) * (R - 2 * D22) *
						Math.pow((R - 2 * D22) * (R - 2 * D22) + a02 * a02, -1.5) -
				0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + a02 * a02, -1.5);
	}

	public static double QpipiQpipideriv(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau) {
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.125 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) +
				0.125 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) -
				0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5) +
				0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
	}

	public static double QxxQyyderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, -1.5) -
				0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5) +
				0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
	}

	public static double QpipiQzzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau) {
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.125 / R * (xB[tau] - xA[tau]) * (R - 2 * D22) *
				Math.pow((R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.125 / R * (xB[tau] - xA[tau]) * (R + 2 * D22) *
						Math.pow((R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R - 2 * D22) *
						Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R + 2 * D22) *
						Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5) -
				0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
	}

	public static double QzzQzzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.0625 / R * (xB[tau] - xA[tau]) * (R + 2 * D21 - 2 * D22) *
				Math.pow((R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) + a22 * a22, -1.5) +
				0.0625 / R * (xB[tau] - xA[tau]) * (R + 2 * D21 + 2 * D22) *
						Math.pow((R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) + a22 * a22, -1.5) +
				0.0625 / R * (xB[tau] - xA[tau]) * (R - 2 * D21 - 2 * D22) *
						Math.pow((R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) + a22 * a22, -1.5) +
				0.0625 / R * (xB[tau] - xA[tau]) * (R - 2 * D21 + 2 * D22) *
						Math.pow((R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R + 2 * D21) *
						Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R - 2 * D21) *
						Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R + 2 * D22) *
						Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R - 2 * D22) *
						Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5) +
				0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
	}

	public static double QpizQpizderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau) {
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.125 / R * (xB[tau] - xA[tau]) * (R + D21 - D22) *
				Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R + D21 - D22) *
						Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R + D21 + D22) *
						Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) +
				0.125 / R * (xB[tau] - xA[tau]) * (R + D21 + D22) *
						Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R - D21 - D22) *
						Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) +
				0.125 / R * (xB[tau] - xA[tau]) * (R - D21 - D22) *
						Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) +
				0.125 / R * (xB[tau] - xA[tau]) * (R - D21 + D22) *
						Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R - D21 + D22) *
						Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5);
	}

	public static double ssssderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double ssppippideriv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double sspzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double ppippissderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}

	public static double pzpzssderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}

	public static double ppippippippideriv(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										   int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) +
				QpipiQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double pxpxpypyderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) +
				QxxQyyderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double ppippipzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) +
				QpipiQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double pzpzppippideriv(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) +
				QpipiQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}

	public static double pzpzpzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau) {
		return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) +
				qQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) +
				QzzQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double spzssderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return -quzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}

	public static double spzppippideriv(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau) {
		return -quzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) +
				uzQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double spzpzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									  int tau) {
		return -quzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) +
				uzQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double ssspzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return quzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double ppippispzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau) {
		return quzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) -
				uzQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}

	public static double pzpzspzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									  int tau) {
		return quzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) -
				uzQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}

	public static double sppisppideriv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau) {
		return upiupideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double spzspzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return uzuzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double sppippipzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau) {
		return upiQpizderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double ppipzsppideriv(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau) {
		return -upiQpizderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}

	public static double ppipzppipzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau) {
		return QpizQpizderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}

	public static double pxpypxpyderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau) {
		return 0.5 * (ppippippippideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) -
				pxpxpypyderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau));
	}

	public static double upiupideriv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		double R = GTO.R(xA, xB);
		double a11 = p11 + p12;
		return 0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + (D11 - D12) * (D11 - D12) + a11 * a11, -1.5) -
				0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + (D11 + D12) * (D11 + D12) + a11 * a11, -1.5);
	}

	public static double uzuzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		double R = GTO.R(xA, xB);
		double a11 = p11 + p12;
		return 0.25 / R * (xB[tau] - xA[tau]) * (R + D11 - D12) *
				Math.pow((R + D11 - D12) * (R + D11 - D12) + a11 * a11, -1.5) -
				0.25 / R * (xB[tau] - xA[tau]) * (R + D11 + D12) *
						Math.pow((R + D11 + D12) * (R + D11 + D12) + a11 * a11, -1.5) -
				0.25 / R * (xB[tau] - xA[tau]) * (R - D11 - D12) *
						Math.pow((R - D11 - D12) * (R - D11 - D12) + a11 * a11, -1.5) +
				0.25 / R * (xB[tau] - xA[tau]) * (R - D11 + D12) *
						Math.pow((R - D11 + D12) * (R - D11 + D12) + a11 * a11, -1.5);
	}

	public static double upiQpizderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									  int tau) {
		double R = GTO.R(xA, xB);
		double a12 = p11 + p22;
		return -0.25 / R * (xB[tau] - xA[tau]) * (R - D22) *
				Math.pow((R - D22) * (R - D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -1.5) +
				0.25 / R * (xB[tau] - xA[tau]) * (R - D22) *
						Math.pow((R - D22) * (R - D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -1.5) +
				0.25 / R * (xB[tau] - xA[tau]) * (R + D22) *
						Math.pow((R + D22) * (R + D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -1.5) -
				0.25 / R * (xB[tau] - xA[tau]) * (R + D22) *
						Math.pow((R + D22) * (R + D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -1.5);
	}

	public static double uzQpipideriv(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									  int tau) {
		double R = GTO.R(xA, xB);
		double a12 = p11 + p22;
		return -0.25 / R * (xB[tau] - xA[tau]) * (R + D11) *
				Math.pow((R + D11) * (R + D11) + 4 * D22 * D22 + a12 * a12, -1.5) +
				0.25 / R * (xB[tau] - xA[tau]) * (R - D11) *
						Math.pow((R - D11) * (R - D11) + 4 * D22 * D22 + a12 * a12, -1.5) +
				0.25 / R * (xB[tau] - xA[tau]) * (R + D11) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5) -
				0.25 / R * (xB[tau] - xA[tau]) * (R - D11) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
	}

	public static double uzQzzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		double R = GTO.R(xA, xB);
		double a12 = p11 + p22;
		return -0.125 / R * (xB[tau] - xA[tau]) * (R + D11 - 2 * D22) *
				Math.pow((R + D11 - 2 * D22) * (R + D11 - 2 * D22) + a12 * a12, -1.5) +
				0.125 / R * (xB[tau] - xA[tau]) * (R - D11 - 2 * D22) *
						Math.pow((R - D11 - 2 * D22) * (R - D11 - 2 * D22) + a12 * a12, -1.5) -
				0.125 / R * (xB[tau] - xA[tau]) * (R + D11 + 2 * D22) *
						Math.pow((R + D11 + 2 * D22) * (R + D11 + 2 * D22) + a12 * a12, -1.5) +
				0.125 / R * (xB[tau] - xA[tau]) * (R - D11 + 2 * D22) *
						Math.pow((R - D11 + 2 * D22) * (R - D11 + 2 * D22) + a12 * a12, -1.5) +
				0.25 / R * (xB[tau] - xA[tau]) * (R + D11) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5) -
				0.25 / R * (xB[tau] - xA[tau]) * (R - D11) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
	}

	public static double qq(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
							double p1deriv, double p2deriv) {
		return 0;
	}

	public static double quz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							 double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
							 double p1deriv, double p2deriv) {
		double a01 = p01 + p12;
		if (num == 0) return 0;
		if (num == 1)
			return -0.5 * ((R + D12) * D1deriv + a01 * p1deriv) * Math.pow((R + D12) * (R + D12) + a01 * a01, -1.5) +
					0.5 * ((D12 - R) * D1deriv + a01 * p1deriv) * Math.pow((R - D12) * (R - D12) + a01 * a01, -1.5);
		if (num == 2)
			return -0.5 * ((R + D11) * D1deriv + a01 * p1deriv) * Math.pow((R + D11) * (R + D11) + a01 * a01, -1.5) +
					0.5 * ((D11 - R) * D1deriv + a01 * p1deriv) * Math.pow((R - D11) * (R - D11) + a01 * a01, -1.5);
		return 0;
	}

	public static double qQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								double p1deriv, double p2deriv) {
		double a02 = p01 + p22;
		if (num == 0) return 0;
		if (num == 1)
			return -0.5 * (4 * D22 * D2deriv + a02 * p2deriv) * Math.pow(R * R + 4 * D22 * D22 + a02 * a02, -1.5) +
					0.5 * a02 * p2deriv * Math.pow(R * R + a02 * a02, -1.5);
		if (num == 2)
			return -0.5 * (4 * D21 * D2deriv + a02 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a02 * a02, -1.5) +
					0.5 * a02 * p2deriv * Math.pow(R * R + a02 * a02, -1.5);
		return 0;
	}

	public static double qQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							  double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
							  double p1deriv, double p2deriv) {
		double a02 = p01 + p22;
		if (num == 0) return 0;
		if (num == 1) return -0.25 * (2 * (R + 2 * D22) * D2deriv + a02 * p2deriv) *
				Math.pow((R + 2 * D22) * (R + 2 * D22) + a02 * a02, -1.5) -
				0.25 * (2 * (2 * D22 - R) * D2deriv + a02 * p2deriv) *
						Math.pow((R - 2 * D22) * (R - 2 * D22) + a02 * a02, -1.5) +
				0.5 * a02 * p2deriv * Math.pow(R * R + a02 * a02, -1.5);
		if (num == 2) return -0.25 * (2 * (R + 2 * D21) * D2deriv + a02 * p2deriv) *
				Math.pow((R + 2 * D21) * (R + 2 * D21) + a02 * a02, -1.5) -
				0.25 * (2 * (2 * D21 - R) * D2deriv + a02 * p2deriv) *
						Math.pow((R - 2 * D21) * (R - 2 * D21) + a02 * a02, -1.5) +
				0.5 * a02 * p2deriv * Math.pow(R * R + a02 * a02, -1.5);
		return 0;
	}

	public static double upiupi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								double p1deriv, double p2deriv) {
		double a11 = p11 + p12;
		if (num == 0) return -0.5 * ((D11 - D12) * D1deriv + a11 * p1deriv) *
				Math.pow(R * R + (D11 - D12) * (D11 - D12) + a11 * a11, -1.5) +
				0.5 * ((D11 + D12) * D1deriv + a11 * p1deriv) *
						Math.pow(R * R + (D11 + D12) * (D11 + D12) + a11 * a11, -1.5);
		if (num == 1) return -0.5 * ((D12 - D11) * D1deriv + a11 * p1deriv) *
				Math.pow(R * R + (D11 - D12) * (D11 - D12) + a11 * a11, -1.5) +
				0.5 * ((D11 + D12) * D1deriv + a11 * p1deriv) *
						Math.pow(R * R + (D11 + D12) * (D11 + D12) + a11 * a11, -1.5);
		if (num == 2) return -0.5 * (2 * a11 * p1deriv) * Math.pow(R * R + a11 * a11, -1.5) +
				0.5 * (4 * D11 * D1deriv + 2 * a11 * p1deriv) * Math.pow(R * R + 4 * D11 * D11 + a11 * a11, -1.5);
		return 0;
	}

	public static double uzuz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							  double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
							  double p1deriv, double p2deriv) {
		double a11 = p11 + p12;
		if (num == 0) return -0.25 * ((R + D11 - D12) * D1deriv + a11 * p1deriv) *
				Math.pow((R + D11 - D12) * (R + D11 - D12) + a11 * a11, -1.5) +
				0.25 * ((R + D11 + D12) * D1deriv + a11 * p1deriv) *
						Math.pow((R + D11 + D12) * (R + D11 + D12) + a11 * a11, -1.5) +
				0.25 * ((D11 + D12 - R) * D1deriv + a11 * p1deriv) *
						Math.pow((R - D11 - D12) * (R - D11 - D12) + a11 * a11, -1.5) -
				0.25 * ((D11 - D12 - R) * D1deriv + a11 * p1deriv) *
						Math.pow((R - D11 + D12) * (R - D11 + D12) + a11 * a11, -1.5);
		if (num == 1) return -0.25 * ((D12 - D11 - R) * D1deriv + a11 * p1deriv) *
				Math.pow((R + D11 - D12) * (R + D11 - D12) + a11 * a11, -1.5) +
				0.25 * ((R + D11 + D12) * D1deriv + a11 * p1deriv) *
						Math.pow((R + D11 + D12) * (R + D11 + D12) + a11 * a11, -1.5) +
				0.25 * ((D11 + D12 - R) * D1deriv + a11 * p1deriv) *
						Math.pow((R - D11 - D12) * (R - D11 - D12) + a11 * a11, -1.5) -
				0.25 * ((D12 + R - D11) * D1deriv + a11 * p1deriv) *
						Math.pow((R - D11 + D12) * (R - D11 + D12) + a11 * a11, -1.5);
		if (num == 2) return -a11 * p1deriv * Math.pow(R * R + a11 * a11, -1.5) +
				0.25 * (2 * (R + 2 * D11) * D1deriv + 2 * a11 * p1deriv) *
						Math.pow((R + 2 * D11) * (R + 2 * D11) + a11 * a11, -1.5) +
				0.25 * (2 * (2 * D11 - R) * D1deriv + 2 * a11 * p1deriv) *
						Math.pow((R - 2 * D11) * (R - 2 * D11) + a11 * a11, -1.5);
		return 0;
	}

	public static double upiQpiz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								 double p1deriv, double p2deriv) {
		double a12 = p11 + p22;
		if (num == 0) return 0.25 * ((D11 - D22) * D1deriv + a12 * p1deriv) *
				Math.pow((R - D22) * (R - D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -1.5) -
				0.25 * ((D11 + D22) * D1deriv + a12 * p1deriv) *
						Math.pow((R - D22) * (R - D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -1.5) -
				0.25 * ((D11 - D22) * D1deriv + a12 * p1deriv) *
						Math.pow((R + D22) * (R + D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -1.5) +
				0.25 * ((D11 + D22) * D1deriv + a12 * p1deriv) *
						Math.pow((R + D22) * (R + D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -1.5);
		if (num == 1) return 0.25 * ((2 * D22 - R - D11) * D2deriv + a12 * p2deriv) *
				Math.pow((R - D22) * (R - D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -1.5) -
				0.25 * ((2 * D22 + D11 - R) * D2deriv + a12 * p2deriv) *
						Math.pow((R - D22) * (R - D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -1.5) -
				0.25 * ((2 * D22 + R - D11) * D2deriv + a12 * p2deriv) *
						Math.pow((R + D22) * (R + D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -1.5) +
				0.25 * ((2 * D22 + D11 + R) * D2deriv + a12 * p2deriv) *
						Math.pow((R + D22) * (R + D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -1.5);
		if (num == 2)
			return 0.25 * ((D21 - R) * D2deriv + (D11 - D21) * (D1deriv - D2deriv) + a12 * (p1deriv + p2deriv)) *
					Math.pow((R - D21) * (R - D21) + (D11 - D21) * (D11 - D21) + a12 * a12, -1.5) -
					0.25 * ((D21 - R) * D2deriv + (D11 + D21) * (D1deriv + D2deriv) + a12 * (p1deriv + p2deriv)) *
							Math.pow((R - D21) * (R - D21) + (D11 + D21) * (D11 + D21) + a12 * a12, -1.5) -
					0.25 * ((D21 + R) * D2deriv + (D11 - D21) * (D1deriv - D2deriv) + a12 * (p1deriv + p2deriv)) *
							Math.pow((R + D21) * (R + D21) + (D11 - D21) * (D11 - D21) + a12 * a12, -1.5) +
					0.25 * ((D21 + R) * D2deriv + (D11 + D21) * (D1deriv + D2deriv) + a12 * (p1deriv + p2deriv)) *
							Math.pow((R + D21) * (R + D21) + (D11 + D21) * (D11 + D21) + a12 * a12, -1.5);
		return 0;
	}

	public static double uzQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								 double p1deriv, double p2deriv) {
		double a12 = p11 + p22;
		if (num == 0) return 0.25 * ((R + D11) * D1deriv + a12 * p1deriv) *
				Math.pow((R + D11) * (R + D11) + 4 * D22 * D22 + a12 * a12, -1.5) -
				0.25 * ((D11 - R) * D1deriv + a12 * p1deriv) *
						Math.pow((R - D11) * (R - D11) + 4 * D22 * D22 + a12 * a12, -1.5) -
				0.25 * ((R + D11) * D1deriv + a12 * p1deriv) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5) +
				0.25 * ((D11 - R) * D1deriv + a12 * p1deriv) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		if (num == 1) return 0.25 * (4 * D22 * D2deriv + a12 * p2deriv) *
				Math.pow((R + D11) * (R + D11) + 4 * D22 * D22 + a12 * a12, -1.5) -
				0.25 * (4 * D22 * D2deriv + a12 * p2deriv) *
						Math.pow((R - D11) * (R - D11) + 4 * D22 * D22 + a12 * a12, -1.5) -
				0.25 * (a12 * p2deriv) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5) +
				0.25 * (a12 * p2deriv) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		if (num == 2) return 0.25 * ((R + D11) * D1deriv + 4 * D21 * D2deriv + a12 * (p1deriv + p2deriv)) *
				Math.pow((R + D11) * (R + D11) + 4 * D21 * D21 + a12 * a12, -1.5) -
				0.25 * ((D11 - R) * D1deriv + 4 * D21 * D2deriv + a12 * (p1deriv + p2deriv)) *
						Math.pow((R - D11) * (R - D11) + 4 * D21 * D21 + a12 * a12, -1.5) -
				0.25 * ((R + D11) * D1deriv + a12 * (p1deriv + p2deriv)) *
						Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5) +
				0.25 * ((D11 - R) * D1deriv + a12 * (p1deriv + p2deriv)) *
						Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		return 0;
	}

	public static double uzQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							   double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
							   double p1deriv, double p2deriv) {
		double a12 = p11 + p22;
		if (num == 0) return +0.125 * ((R + D11 - 2 * D22) * D1deriv + a12 * p1deriv) *
				Math.pow((R + D11 - 2 * D22) * (R + D11 - 2 * D22) + a12 * a12, -1.5) -
				0.125 * ((D11 + 2 * D22 - R) * D1deriv + a12 * p1deriv) *
						Math.pow((R - D11 - 2 * D22) * (R - D11 - 2 * D22) + a12 * a12, -1.5) +
				0.125 * ((R + D11 + 2 * D22) * D1deriv + a12 * p1deriv) *
						Math.pow((R + D11 + 2 * D22) * (R + D11 + 2 * D22) + a12 * a12, -1.5) -
				0.125 * ((D11 - 2 * D22 - R) * D1deriv + a12 * p1deriv) *
						Math.pow((R - D11 + 2 * D22) * (R - D11 + 2 * D22) + a12 * a12, -1.5) -
				0.25 * ((R + D11) * D1deriv + a12 * p1deriv) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5) +
				0.25 * ((D11 - R) * D1deriv + a12 * p1deriv) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		if (num == 1) return +0.125 * (2 * (2 * D22 - D11 - R) * D2deriv + a12 * p2deriv) *
				Math.pow((R + D11 - 2 * D22) * (R + D11 - 2 * D22) + a12 * a12, -1.5) -
				0.125 * (2 * (D11 + 2 * D22 - R) * D2deriv + a12 * p2deriv) *
						Math.pow((R - D11 - 2 * D22) * (R - D11 - 2 * D22) + a12 * a12, -1.5) +
				0.125 * (2 * (R + D11 + 2 * D22) * D2deriv + a12 * p2deriv) *
						Math.pow((R + D11 + 2 * D22) * (R + D11 + 2 * D22) + a12 * a12, -1.5) -
				0.125 * (2 * (2 * D22 + R - D11) * D2deriv + a12 * p2deriv) *
						Math.pow((R - D11 + 2 * D22) * (R - D11 + 2 * D22) + a12 * a12, -1.5) -
				0.25 * (a12 * p2deriv) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5) +
				0.25 * (a12 * p2deriv) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		if (num == 2) return +0.125 * ((R + D11 - 2 * D21) * (D1deriv - 2 * D2deriv) + a12 * (p1deriv + p2deriv)) *
				Math.pow((R + D11 - 2 * D21) * (R + D11 - 2 * D21) + a12 * a12, -1.5) -
				0.125 * ((D11 + 2 * D21 - R) * (D1deriv + 2 * D2deriv) + a12 * (p1deriv + p2deriv)) *
						Math.pow((R - D11 - 2 * D21) * (R - D11 - 2 * D21) + a12 * a12, -1.5) +
				0.125 * ((R + D11 + 2 * D21) * (D1deriv + 2 * D2deriv) + a12 * (p1deriv + p2deriv)) *
						Math.pow((R + D11 + 2 * D21) * (R + D11 + 2 * D21) + a12 * a12, -1.5) -
				0.125 * ((R - D11 + 2 * D21) * (2 * D2deriv - D1deriv) + a12 * (p1deriv + p2deriv)) *
						Math.pow((R - D11 + 2 * D21) * (R - D11 + 2 * D21) + a12 * a12, -1.5) -
				0.25 * ((R + D11) * D1deriv + a12 * (p1deriv + p2deriv)) *
						Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5) +
				0.25 * ((D11 - R) * D1deriv + a12 * (p1deriv + p2deriv)) *
						Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		return 0;
	}

	public static double QpipiQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) return -0.125 * (4 * (D21 - D22) * D2deriv + a22 * p2deriv) *
				Math.pow(R * R + 4 * (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) -
				0.125 * (4 * (D21 + D22) * D2deriv + a22 * p2deriv) *
						Math.pow(R * R + 4 * (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) +
				0.25 * (4 * D21 * D2deriv + a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.25 * (a22 * p2deriv) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5) -
				0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		if (num == 1) return -0.125 * (4 * (D22 - D21) * D2deriv + a22 * p2deriv) *
				Math.pow(R * R + 4 * (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) -
				0.125 * (4 * (D21 + D22) * D2deriv + a22 * p2deriv) *
						Math.pow(R * R + 4 * (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) +
				0.25 * (a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.25 * (4 * D22 * D2deriv + a22 * p2deriv) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5) -
				0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		if (num == 2) return -0.75 * a22 * p2deriv * Math.pow(R * R + a22 * a22, -1.5) -
				0.125 * (16 * D21 * D2deriv + 2 * a22 * p2deriv) * Math.pow(R * R + 16 * D21 * D21 + a22 * a22, -1.5) +
				0.5 * (4 * D21 * D2deriv + 2 * a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5);
		return 0;
	}

	public static double QxxQyy(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								double p1deriv, double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) return -0.25 * (4 * D21 * D2deriv + a22 * p2deriv) *
				Math.pow(R * R + 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, -1.5) +
				0.25 * (4 * D21 * D2deriv + a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.25 * (a22 * p2deriv) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5) -
				0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		if (num == 1) return -0.25 * (4 * D22 * D2deriv + a22 * p2deriv) *
				Math.pow(R * R + 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, -1.5) +
				0.25 * (a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.25 * (4 * D22 * D2deriv + a22 * p2deriv) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5) -
				0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		if (num == 2)
			return -0.25 * (8 * D22 * D2deriv + 2 * a22 * p2deriv) * Math.pow(R * R + 8 * D21 * D21 + a22 * a22,
					-1.5) +
					0.5 * (4 * D22 * D2deriv + 2 * a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) -
					0.25 * (2 * a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		return 0;
	}

	public static double QpipiQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, int num, double D1deriv,
								  double D2deriv,
								  double p1deriv, double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) return -0.125 * (4 * D21 * D2deriv + a22 * p2deriv) *
				Math.pow((R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.125 * (4 * D21 * D2deriv + a22 * p2deriv) *
						Math.pow((R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.125 * (a22 * p2deriv) * Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5) +
				0.125 * (a22 * p2deriv) * Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5) +
				0.25 * (4 * D21 * D2deriv + a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		if (num == 1) return -0.125 * (2 * (2 * D22 - R) * D2deriv + a22 * p2deriv) *
				Math.pow((R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.125 * (2 * (2 * D22 + R) * D2deriv + a22 * p2deriv) *
						Math.pow((R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.125 * (2 * (2 * D22 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5) +
				0.125 * (2 * (2 * D22 + R) * D2deriv + a22 * p2deriv) *
						Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5) +
				0.25 * (a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		if (num == 2) return -0.125 * ((2 * (2 * D21 - R) + 4 * D21) * D2deriv + 2 * a22 * p2deriv) *
				Math.pow((R - 2 * D21) * (R - 2 * D21) + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.125 * ((2 * (2 * D21 + R) + 4 * D21) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R + 2 * D21) * (R + 2 * D21) + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.125 * (2 * (2 * D21 - R) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5) +
				0.125 * (2 * (2 * D21 + R) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5) +
				0.25 * (4 * D21 * D2deriv + 2 * a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.25 * (2 * a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		return 0;
	}

	public static double QzzQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								double p1deriv, double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) return -0.0625 * (2 * (R + 2 * D21 - 2 * D22) * D2deriv + a22 * p2deriv) *
				Math.pow((R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) + a22 * a22, -1.5) -
				0.0625 * (2 * (R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
						Math.pow((R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) + a22 * a22, -1.5) -
				0.0625 * (2 * (2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) + a22 * a22, -1.5) -
				0.0625 * (2 * (2 * D21 - R - 2 * D22) * D2deriv + a22 * p2deriv) *
						Math.pow((R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) + a22 * a22, -1.5) +
				0.125 * (2 * (R + 2 * D21) * D2deriv + a22 * p2deriv) *
						Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5) +
				0.125 * (2 * (2 * D21 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5) +
				0.125 * (a22 * p2deriv) * Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5) +
				0.125 * (a22 * p2deriv) * Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5) -
				0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		if (num == 1) return -0.0625 * (2 * (2 * D22 - R - 2 * D21) * D2deriv + a22 * p2deriv) *
				Math.pow((R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) + a22 * a22, -1.5) -
				0.0625 * (2 * (R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
						Math.pow((R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) + a22 * a22, -1.5) -
				0.0625 * (2 * (2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) + a22 * a22, -1.5) -
				0.0625 * (2 * (R - 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
						Math.pow((R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) + a22 * a22, -1.5) +
				0.125 * (a22 * p2deriv) * Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5) +
				0.125 * (a22 * p2deriv) * Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5) +
				0.125 * (2 * (R + 2 * D22) * D2deriv + a22 * p2deriv) *
						Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5) +
				0.125 * (2 * (2 * D22 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5) -
				0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		if (num == 2) return -0.0625 * (4 * (R + 4 * D21) * D2deriv + 2 * a22 * p2deriv) *
				Math.pow((R + 4 * D21) * (R + 4 * D21) + a22 * a22, -1.5) -
				0.0625 * (4 * (4 * D21 - R) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R - 4 * D21) * (R - 4 * D21) + a22 * a22, -1.5) +
				0.25 * (2 * (R + 2 * D21) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5) +
				0.25 * (2 * (2 * D21 - R) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5) -
				0.75 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		return 0;
	}

	public static double QpizQpiz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, int num, double D1deriv,
								  double D2deriv,
								  double p1deriv, double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) return -0.125 * ((R + 2 * D21 - 2 * D22) * D2deriv + a22 * p2deriv) *
				Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) +
				0.125 * ((R + 2 * D21) * D2deriv + a22 * p2deriv) *
						Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) +
				0.125 * ((R + 2 * D21) * D2deriv + a22 * p2deriv) *
						Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) -
				0.125 * ((R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
						Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) +
				0.125 * ((2 * D21 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) -
				0.125 * ((2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) -
				0.125 * ((2 * D21 - 2 * D22 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) +
				0.125 * ((2 * D21 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5);
		if (num == 1) return -0.125 * ((2 * D22 - 2 * D21 - R) * D2deriv + a22 * p2deriv) *
				Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) +
				0.125 * ((2 * D22 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) +
				0.125 * ((R + 2 * D22) * D2deriv + a22 * p2deriv) *
						Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) -
				0.125 * ((R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
						Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) +
				0.125 * ((2 * D22 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) -
				0.125 * ((2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
						Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5) -
				0.125 * ((2 * D22 + R - 2 * D21) * D2deriv + a22 * p2deriv) *
						Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5) +
				0.125 * ((R + 2 * D22) * D2deriv + a22 * p2deriv) *
						Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5);
		if (num == 2) return -0.5 * a22 * p2deriv * Math.pow(R * R + a22 * a22, -1.5) +
				0.25 * (4 * D21 * D2deriv + 2 * a22 * p2deriv) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5) +
				0.125 * (2 * (R + 2 * D21) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5) +
				0.125 * (2 * (2 * D21 - R) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5) -
				0.125 * ((2 * (2 * D21 + R) + 4 * D21) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R + 2 * D21) * (R + 2 * D21) + 4 * D21 * D21 + a22 * a22, -1.5) -
				0.125 * ((2 * (2 * D21 - R) + 4 * D21) * D2deriv + 2 * a22 * p2deriv) *
						Math.pow((R - 2 * D21) * (R - 2 * D21) + 4 * D21 * D21 + a22 * a22, -1.5);
		return 0;
	}

	public static double ssssderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, int num, double D1deriv,
								   double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ssppippideriv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double sspzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippissderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzssderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippippippideriv(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, int num,
										   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv) +
				QpipiQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pxpxpypyderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv) +
				QxxQyy(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippipzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv) +
				QpipiQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double pzpzppippideriv(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				QpipiQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzpzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				QzzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzssderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzppippideriv(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				uzQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzpzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, int num,
									  double D1deriv,
									  double D2deriv, double p1deriv, double p2deriv) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				uzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ssspzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippispzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) -
				uzQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzspzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, int num,
									  double D1deriv,
									  double D2deriv, double p1deriv, double p2deriv) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) -
				uzQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double sppisppideriv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return upiupi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzspzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return uzuz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double sppippipzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return upiQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppipzsppideriv(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return -upiQpiz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
				p2deriv);
	}

	public static double ppipzppipzderiv(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return QpizQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double pxpypxpyderiv(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return 0.5 *
				(ppippippippideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv) -
						pxpxpypyderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv,
								p1deriv, p2deriv));
	}

	public static int f(int num) {
		if (num == 0 || num == 1) return 1 - num;
		else if (num == 2) return num;
		return -1;
	}

	private static double qqderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
		double sum = 0;
		double a00 = p01 + p02;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum = generalizedform2(0, a00 * a00, R);
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a00 * a00, R);
	}

	private static double quzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
		double sum = 0;
		double a01 = p01 + p12;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum =
				generalizedform2(D12, a01 * a01, R) * 0.5 + generalizedform2(-D12, a01 * a01, R) * -0.5;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(D12, a01 * a01, R) * 0.5 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D12, a01 * a01, R) * -0.5;
	}

	private static double qQpipideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		double sum = 0;
		double a02 = p01 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2)
			sum = generalizedform2(0, 4 * D22 * D22 + a02 * a02, R) * 0.5 + generalizedform2(0, a02 * a02, R) * -0.5;
		return sum +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D22 * D22 + a02 * a02, R) * 0.5 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a02 * a02, R) * -0.5;
	}

	private static double qQzzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau1,
									 int tau2) {
		double sum = 0;
		double a02 = p01 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum =
				generalizedform2(2 * D22, a02 * a02, R) * 0.25 + generalizedform2(-2 * D22, a02 * a02, R) * 0.25 +
						generalizedform2(0, a02 * a02, R) * -0.5;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(2 * D22, a02 * a02, R) * 0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D22, a02 * a02, R) * 0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a02 * a02, R) * -0.5;
	}

	private static double upiupideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		double sum = 0;
		double a11 = p11 + p12;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum = generalizedform2(0, (D11 - D12) * (D11 - D12) + a11 * a11, R) * 0.5 +
				generalizedform2(0, (D11 + D12) * (D11 + D12) + a11 * a11, R) * -0.5;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, (D11 - D12) * (D11 - D12) + a11 * a11, R) * 0.5 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
						generalizedform(0, (D11 + D12) * (D11 + D12) + a11 * a11, R) * -0.5;
	}

	private static double uzuzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau1,
									 int tau2) {
		double sum = 0;
		double a11 = p11 + p12;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum =
				generalizedform2(D11 - D12, a11 * a11, R) * 0.25 + generalizedform2(D11 + D12, a11 * a11, R) * -0.25 +
						generalizedform2(-D11 - D12, a11 * a11, R) * -0.25 +
						generalizedform2(-D11 + D12, a11 * a11, R) * 0.25;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(D11 - D12, a11 * a11, R) * 0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(D11 + D12, a11 * a11, R) * -0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11 - D12, a11 * a11, R) * -0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11 + D12, a11 * a11, R) * 0.25;
	}

	private static double upiQpizderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		double sum = 0;
		double a12 = p11 + p22;
		double R = GTO.R(xA, xB);
		double denom = (D11 - D22) * (D11 - D22) + a12 * a12;
		double denom2 = (D11 + D22) * (D11 + D22) + a12 * a12;
		if (tau1 == tau2) sum = generalizedform2(-D22, denom, R) * -0.25 + generalizedform2(-D22, denom2, R) * 0.25 +
				generalizedform2(+D22, denom, R) * 0.25 + generalizedform2(+D22, denom2, R) * -0.25;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D22, denom, R) * -0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D22, denom2, R) * 0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D22, denom, R) * 0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D22, denom2, R) * -0.25;
	}

	private static double uzQpipideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		double sum = 0;
		double a12 = p11 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum = generalizedform2(+D11, 4 * D22 * D22 + a12 * a12, R) * -0.25 +
				generalizedform2(-D11, 4 * D22 * D22 + a12 * a12, R) * 0.25 +
				generalizedform2(+D11, a12 * a12, R) * 0.25 + generalizedform2(-D11, a12 * a12, R) * -0.25;
		return sum +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11, 4 * D22 * D22 + a12 * a12, R) *
						-0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11, 4 * D22 * D22 + a12 * a12, R) *
						0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11, a12 * a12, R) * 0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11, a12 * a12, R) * -0.25;
	}

	private static double uzQzzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									  int tau1, int tau2) {
		double sum = 0;
		double a12 = p11 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum = generalizedform2(+D11 - 2 * D22, a12 * a12, R) * -0.125 +
				generalizedform2(-D11 - 2 * D22, a12 * a12, R) * 0.125 +
				generalizedform2(+D11 + 2 * D22, a12 * a12, R) * -0.125 +
				generalizedform2(-D11 + 2 * D22, a12 * a12, R) * 0.125 + generalizedform2(+D11, a12 * a12, R) * 0.25 +
				generalizedform2(-D11, a12 * a12, R) * -0.25;
		return sum +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11 - 2 * D22, a12 * a12, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11 - 2 * D22, a12 * a12, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11 + 2 * D22, a12 * a12, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11 + 2 * D22, a12 * a12, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11, a12 * a12, R) * 0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11, a12 * a12, R) * -0.25;
	}

	private static double QpipiQpipideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										   int tau1, int tau2) {
		double sum = 0;
		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum = generalizedform2(0, 4 * (D21 - D22) * (D21 - D22) + a22 * a22, R) * 0.125 +
				generalizedform2(0, 4 * (D21 + D22) * (D21 + D22) + a22 * a22, R) * 0.125 +
				generalizedform2(0, 4 * D21 * D21 + a22 * a22, R) * -0.25 +
				generalizedform2(0, 4 * D22 * D22 + a22 * a22, R) * -0.25 + generalizedform2(0, a22 * a22, R) * 0.25;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * (D21 - D22) * (D21 - D22) + a22 * a22, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
						generalizedform(0, 4 * (D21 + D22) * (D21 + D22) + a22 * a22, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D21 * D21 + a22 * a22, R) *
						-0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D22 * D22 + a22 * a22, R) *
						-0.25 + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a22 * a22, R) * 0.25;
	}

	private static double QxxQyyderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		double sum = 0;
		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum = generalizedform2(0, 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, R) * 0.25 +
				generalizedform2(0, 4 * D21 * D21 + a22 * a22, R) * -0.25 +
				generalizedform2(0, 4 * D22 * D22 + a22 * a22, R) * -0.25 + generalizedform2(0, a22 * a22, R) * 0.25;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, R) * 0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D21 * D21 + a22 * a22, R) *
						-0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D22 * D22 + a22 * a22, R) *
						-0.25 + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a22 * a22, R) * 0.25;
	}

	private static double QpipiQzzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		double sum = 0;
		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);
		double denom = 4 * D21 * D21 + a22 * a22;
		if (tau1 == tau2) sum =
				generalizedform2(-2 * D22, denom, R) * 0.125 + generalizedform2(+2 * D22, denom, R) * 0.125 +
						generalizedform2(-2 * D22, a22 * a22, R) * -0.125 +
						generalizedform2(+2 * D22, a22 * a22, R) * -0.125 + generalizedform2(0, denom, R) * -0.25 +
						generalizedform2(0, a22 * a22, R) * 0.25;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D22, denom, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D22, denom, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D22, a22 * a22, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D22, a22 * a22, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, denom, R) * -0.25 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a22 * a22, R) * 0.25;
	}

	private static double QzzQzzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		double sum = 0;
		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) sum = generalizedform2(+2 * D21 - 2 * D22, a22 * a22, R) * 0.0625 +
				generalizedform2(+2 * D21 + 2 * D22, a22 * a22, R) * 0.0625 +
				generalizedform2(-2 * D21 - 2 * D22, a22 * a22, R) * 0.0625 +
				generalizedform2(-2 * D21 + 2 * D22, a22 * a22, R) * 0.0625 +
				generalizedform2(+2 * D21, a22 * a22, R) * -0.125 + generalizedform2(-2 * D21, a22 * a22, R) * -0.125 +
				generalizedform2(+2 * D22, a22 * a22, R) * -0.125 + generalizedform2(-2 * D22, a22 * a22, R) * -0.125 +
				generalizedform2(0, a22 * a22, R) * 0.25;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D21 - 2 * D22, a22 * a22,
				R) *
				0.0625 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D21 + 2 * D22, a22 * a22, R) *
						0.0625 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D21 - 2 * D22, a22 * a22, R) *
						0.0625 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D21 + 2 * D22, a22 * a22, R) *
						0.0625 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D21, a22 * a22, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D21, a22 * a22, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D22, a22 * a22, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D22, a22 * a22, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a22 * a22, R) * 0.25;
	}

	private static double QpizQpizderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		double sum = 0;
		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);
		double denom1 = (D21 - D22) * (D21 - D22) + a22 * a22;
		double denom2 = (D21 + D22) * (D21 + D22) + a22 * a22;
		if (tau1 == tau2) sum =
				generalizedform2(+D21 - D22, denom1, R) * 0.125 + generalizedform2(+D21 - D22, denom2, R) * -0.125 +
						generalizedform2(+D21 + D22, denom1, R) * -0.125 +
						generalizedform2(+D21 + D22, denom2, R) * 0.125 +
						generalizedform2(-D21 - D22, denom1, R) * -0.125 +
						generalizedform2(-D21 - D22, denom2, R) * 0.125 +
						generalizedform2(-D21 + D22, denom1, R) * 0.125 +
						generalizedform2(-D21 + D22, denom2, R) * -0.125;
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D21 - D22, denom1, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D21 - D22, denom2, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D21 + D22, denom1, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D21 + D22, denom2, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D21 - D22, denom1, R) * -0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D21 - D22, denom2, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D21 + D22, denom1, R) * 0.125 +
				(xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D21 + D22, denom2, R) * -0.125;
	}

	public static double ssssderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double ssppippideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double sspzpzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									  int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double ppippissderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
	}

	public static double pzpzssderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									  int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
	}

	public static double ppippippippideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double[] xA, double[] xB,
											int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) +
				QpipiQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double pxpxpypyderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) +
				QxxQyyderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double ppippipzpzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										  int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) +
				QpipiQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double pzpzppippideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										  int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) +
				QpipiQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
	}

	public static double pzpzpzpzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) +
				qQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) +
				QzzQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double spzssderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
		return -quzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
	}

	public static double spzppippideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return -quzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) +
				uzQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double spzpzpzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		return -quzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) +
				uzQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double ssspzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									 double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
		return quzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double ppippispzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return quzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) -
				uzQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
	}

	public static double pzpzspzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		return quzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) -
				uzQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
	}

	public static double sppisppideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		return upiupideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double spzspzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
									  int tau1, int tau2) {
		return uzuzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double sppippipzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return upiQpizderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double ppipzsppideriv2(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return -upiQpizderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
	}

	public static double ppipzppipzderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										  int tau1, int tau2) {
		return QpizQpizderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
	}

	public static double pxpypxpyderiv2(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		return 0.5 * (ppippippippideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) -
				pxpxpypyderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2));
	}

	private static double generalizedform(double a, double b, double R) {
		double denom = (R + a) * (R + a) + b;
		return 1 / (R * R * Math.pow(denom, 1.5)) * (a / R + 3 * (R + a) * (R + a) / (denom));
	}

	private static double generalizedform2(double a, double b, double R) {
		return -(R + a) / (Math.pow((R + a) * (R + a) + b, 1.5) * R);
	}
}
