package nddo.geometry;

import nddo.NDDO6G;
import nddo.NDDOAtom;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import scf.GTO;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;


public class GeometrySecondDerivative {

	private static double generalizedform(double a, double b, double R) {
		double denom = (R + a) * (R + a) + b;
		return 1 / (R * R * Math.pow(denom, 1.5)) *
				(a / R + 3 * (R + a) * (R + a) / (denom));
	}

	private static double generalizedform2(double a, double b, double R) {
		return -(R + a) / (Math.pow((R + a) * (R + a) + b, 1.5) * R);
	}

	private static double qqderiv2(double p01, double p11, double p21,
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double[] xA, double[] xB,
								   int tau1, int tau2) {
		double sum = 0;

		double a00 = p01 + p02;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {

			sum = generalizedform2(0, a00 * a00, R);
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, a00 * a00, R);
	}

	private static double quzderiv2(double p01, double p11, double p21,
									double D11, double D21, double p02,
									double p12, double p22, double D12,
									double D22, double[] xA, double[] xB,
									int tau1, int tau2) {
		double sum = 0;

		double a01 = p01 + p12;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {
			sum = generalizedform2(D12, a01 * a01, R) * 0.5 +
					generalizedform2(-D12, a01 * a01, R) * -0.5;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(D12, a01 * a01, R) * 0.5
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D12, a01 * a01, R) * -0.5;
	}

	private static double qQpipideriv2(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		double sum = 0;

		double a02 = p01 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {
			sum = generalizedform2(0, 4 * D22 * D22 + a02 * a02, R) * 0.5 +
					generalizedform2(0, a02 * a02, R) * -0.5;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * D22 * D22 + a02 * a02, R) * 0.5
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, a02 * a02, R) * -0.5;
	}

	private static double qQzzderiv2(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double[] xA, double[] xB,
									 int tau1, int tau2) {
		double sum = 0;

		double a02 = p01 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {
			sum = generalizedform2(2 * D22, a02 * a02, R) * 0.25 +
					generalizedform2(-2 * D22, a02 * a02, R) * 0.25 +
					generalizedform2(0, a02 * a02, R) * -0.5;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(2 * D22, a02 * a02, R) * 0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-2 * D22, a02 * a02, R) * 0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, a02 * a02, R) * -0.5;
	}

	private static double upiupideriv2(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		double sum = 0;

		double a11 = p11 + p12;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {
			sum = generalizedform2(0, (D11 - D12) * (D11 - D12) + a11 * a11,
					R) * 0.5 +
					generalizedform2(0, (D11 + D12) * (D11 + D12) + a11 * a11,
							R) * -0.5;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, (D11 - D12) * (D11 - D12) + a11 * a11, R) *
				0.5
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, (D11 + D12) * (D11 + D12) + a11 * a11, R) *
				-0.5;
	}

	private static double uzuzderiv2(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double[] xA, double[] xB,
									 int tau1, int tau2) {
		double sum = 0;

		double a11 = p11 + p12;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {
			sum = generalizedform2(D11 - D12, a11 * a11, R) * 0.25 +
					generalizedform2(D11 + D12, a11 * a11, R) * -0.25 +
					generalizedform2(-D11 - D12, a11 * a11, R) * -0.25 +
					generalizedform2(-D11 + D12, a11 * a11, R) * 0.25;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(D11 - D12, a11 * a11, R) * 0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(D11 + D12, a11 * a11, R) * -0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D11 - D12, a11 * a11, R) * -0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D11 + D12, a11 * a11, R) * 0.25;
	}

	private static double upiQpizderiv2(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		double sum = 0;

		double a12 = p11 + p22;
		double R = GTO.R(xA, xB);
		double denom = (D11 - D22) * (D11 - D22) + a12 * a12;
		double denom2 = (D11 + D22) * (D11 + D22) + a12 * a12;
		if (tau1 == tau2) {
			sum = generalizedform2(-D22, denom, R) * -0.25 +
					generalizedform2(-D22, denom2, R) * 0.25 +
					generalizedform2(+D22, denom, R) * 0.25 +
					generalizedform2(+D22, denom2, R) * -0.25;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D22, denom, R) * -0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D22, denom2, R) * 0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D22, denom, R) * 0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D22, denom2, R) * -0.25;
	}

	private static double uzQpipideriv2(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		double sum = 0;

		double a12 = p11 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {
			sum =
					generalizedform2(+D11, 4 * D22 * D22 + a12 * a12, R) *
							-0.25 +
							generalizedform2(-D11, 4 * D22 * D22 + a12 * a12,
									R) *
									0.25 +
							generalizedform2(+D11, a12 * a12, R) * 0.25 +
							generalizedform2(-D11, a12 * a12, R) * -0.25;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D11, 4 * D22 * D22 + a12 * a12, R) * -0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D11, 4 * D22 * D22 + a12 * a12, R) * 0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D11, a12 * a12, R) * 0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D11, a12 * a12, R) * -0.25;
	}

	private static double uzQzzderiv2(double p01, double p11, double p21,
									  double D11, double D21, double p02,
									  double p12, double p22, double D12,
									  double D22, double[] xA, double[] xB,
									  int tau1, int tau2) {
		double sum = 0;

		double a12 = p11 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {
			sum = generalizedform2(+D11 - 2 * D22, a12 * a12, R) * -0.125 +
					generalizedform2(-D11 - 2 * D22, a12 * a12, R) * 0.125 +
					generalizedform2(+D11 + 2 * D22, a12 * a12, R) * -0.125
					+ generalizedform2(-D11 + 2 * D22, a12 * a12, R) * 0.125 +
					generalizedform2(+D11, a12 * a12, R) * 0.25 +
					generalizedform2(-D11, a12 * a12, R) * -0.25;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D11 - 2 * D22, a12 * a12, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D11 - 2 * D22, a12 * a12, R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D11 + 2 * D22, a12 * a12, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D11 + 2 * D22, a12 * a12, R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D11, a12 * a12, R) * 0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D11, a12 * a12, R) * -0.25;
	}

	private static double QpipiQpipideriv2(double p01, double p11, double p21,
										   double D11, double D21, double p02,
										   double p12, double p22, double D12,
										   double D22, double[] xA,
										   double[] xB,
										   int tau1, int tau2) {
		double sum = 0;

		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {
			sum = generalizedform2(0,
					4 * (D21 - D22) * (D21 - D22) + a22 * a22,
					R) * 0.125 + generalizedform2(0,
					4 * (D21 + D22) * (D21 + D22) + a22 * a22, R) * 0.125
					+
					generalizedform2(0, 4 * D21 * D21 + a22 * a22, R) * -0.25 +
					generalizedform2(0, 4 * D22 * D22 + a22 * a22, R) * -0.25 +
					generalizedform2(0, a22 * a22, R) * 0.25;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * (D21 - D22) * (D21 - D22) + a22 * a22,
						R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * (D21 + D22) * (D21 + D22) + a22 * a22,
						R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * D21 * D21 + a22 * a22, R) * -0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * D22 * D22 + a22 * a22, R) * -0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, a22 * a22, R) * 0.25;
	}

	private static double QxxQyyderiv2(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		double sum = 0;

		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);

		if (tau1 == tau2) {
			sum = generalizedform2(0,
					4 * D21 * D21 + 4 * D22 * D22 + a22 * a22,
					R) * 0.25 +
					generalizedform2(0, 4 * D21 * D21 + a22 * a22, R) * -0.25 +
					generalizedform2(0, 4 * D22 * D22 + a22 * a22, R) * -0.25 +
					generalizedform2(0, a22 * a22, R) * 0.25;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22,
						R) * 0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * D21 * D21 + a22 * a22, R) * -0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, 4 * D22 * D22 + a22 * a22, R) * -0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, a22 * a22, R) * 0.25;
	}

	private static double QpipiQzzderiv2(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		double sum = 0;

		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);
		double denom = 4 * D21 * D21 + a22 * a22;
		if (tau1 == tau2) {
			sum = generalizedform2(-2 * D22, denom, R) * 0.125 +
					generalizedform2(+2 * D22, denom, R) * 0.125 +
					generalizedform2(-2 * D22, a22 * a22, R) * -0.125 +
					generalizedform2(+2 * D22, a22 * a22, R) * -0.125 +
					generalizedform2(0, denom, R) * -0.25 +
					generalizedform2(0, a22 * a22, R) * 0.25;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-2 * D22, denom, R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+2 * D22, denom, R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-2 * D22, a22 * a22, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+2 * D22, a22 * a22, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, denom, R) * -0.25
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, a22 * a22, R) * 0.25;
	}


	private static double QzzQzzderiv2(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		double sum = 0;

		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);
		if (tau1 == tau2) {
			sum = generalizedform2(+2 * D21 - 2 * D22, a22 * a22, R) * 0.0625 +
					generalizedform2(+2 * D21 + 2 * D22, a22 * a22, R) *
							0.0625 +
					generalizedform2(-2 * D21 - 2 * D22, a22 * a22, R) *
							0.0625 +
					generalizedform2(-2 * D21 + 2 * D22, a22 * a22, R) * 0.0625
					+ generalizedform2(+2 * D21, a22 * a22, R) * -0.125 +
					generalizedform2(-2 * D21, a22 * a22, R) * -0.125 +
					generalizedform2(+2 * D22, a22 * a22, R) * -0.125 +
					generalizedform2(-2 * D22, a22 * a22, R) * -0.125 +
					generalizedform2(0, a22 * a22, R) * 0.25;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+2 * D21 - 2 * D22, a22 * a22, R) * 0.0625
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+2 * D21 + 2 * D22, a22 * a22, R) * 0.0625
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-2 * D21 - 2 * D22, a22 * a22, R) * 0.0625
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-2 * D21 + 2 * D22, a22 * a22, R) * 0.0625
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+2 * D21, a22 * a22, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-2 * D21, a22 * a22, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+2 * D22, a22 * a22, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-2 * D22, a22 * a22, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(0, a22 * a22, R) * 0.25;
	}

	private static double QpizQpizderiv2(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		double sum = 0;

		double a22 = p21 + p22;
		double R = GTO.R(xA, xB);
		double denom1 = (D21 - D22) * (D21 - D22) + a22 * a22;
		double denom2 = (D21 + D22) * (D21 + D22) + a22 * a22;
		if (tau1 == tau2) {
			sum = generalizedform2(+D21 - D22, denom1, R) * 0.125 +
					generalizedform2(+D21 - D22, denom2, R) * -0.125 +
					generalizedform2(+D21 + D22, denom1, R) * -0.125 +
					generalizedform2(+D21 + D22, denom2, R) * 0.125
					+ generalizedform2(-D21 - D22, denom1, R) * -0.125 +
					generalizedform2(-D21 - D22, denom2, R) * 0.125 +
					generalizedform2(-D21 + D22, denom1, R) * 0.125 +
					generalizedform2(-D21 + D22, denom2, R) * -0.125;
		}
		return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D21 - D22, denom1, R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D21 - D22, denom2, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D21 + D22, denom1, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(+D21 + D22, denom2, R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D21 - D22, denom1, R) * -0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D21 - D22, denom2, R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D21 + D22, denom1, R) * 0.125
				+ (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) *
				generalizedform(-D21 + D22, denom2, R) * -0.125;
	}

	private static double ssssderiv2(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double[] xA, double[] xB,
									 int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2);
	}

	private static double ssppippideriv2(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) +
				qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2);
	}

	private static double sspzpzderiv2(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) +
				qQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2);
	}

	private static double ppippissderiv2(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) +
				qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
						xA, xB, tau1, tau2);
	}

	private static double pzpzssderiv2(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) +
				qQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
						xA, xB, tau1, tau2);
	}

	private static double ppippippippideriv2(double p01, double p11,
											 double p21,
											 double D11, double D21,
											 double p02,
											 double p12, double p22,
											 double D12,
											 double D22, double[] xA,
											 double[] xB, int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) +
				qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2) +
				qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
						xA, xB, tau1, tau2) +
				QpipiQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12,
						D22, xA, xB, tau1, tau2);
	}

	private static double pxpxpypyderiv2(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) +
				qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2) +
				qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
						xA, xB, tau1, tau2) +
				QxxQyyderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2);
	}

	private static double ppippipzpzderiv2(double p01, double p11, double p21,
										   double D11, double D21, double p02,
										   double p12, double p22, double D12,
										   double D22, double[] xA,
										   double[] xB,
										   int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) +
				qQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2) +
				qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
						xA, xB, tau1, tau2) +
				QpipiQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12,
						D22, xA, xB, tau1, tau2);
	}

	private static double pzpzppippideriv2(double p01, double p11, double p21,
										   double D11, double D21, double p02,
										   double p12, double p22, double D12,
										   double D22, double[] xA,
										   double[] xB,
										   int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) +
				qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2) +
				qQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
						xA, xB, tau1, tau2) +
				QpipiQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11,
						D21, xA, xB, tau1, tau2);
	}

	private static double pzpzpzpzderiv2(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) +
				qQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2) +
				qQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
						xA, xB, tau1, tau2) +
				QzzQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2);
	}

	private static double spzssderiv2(double p01, double p11, double p21,
									  double D11, double D21, double p02,
									  double p12, double p22, double D12,
									  double D22, double[] xA, double[] xB,
									  int tau1, int tau2) {
		return -quzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA,
				xB, tau1, tau2);
	}

	private static double spzppippideriv2(double p01, double p11, double p21,
										  double D11, double D21, double p02,
										  double p12, double p22, double D12,
										  double D22, double[] xA, double[] xB,
										  int tau1, int tau2) {
		return -quzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA,
				xB, tau1, tau2) +
				uzQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2);
	}

	private static double spzpzpzderiv2(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		return -quzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA,
				xB, tau1, tau2) +
				uzQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						xA, xB, tau1, tau2);
	}

	private static double ssspzderiv2(double p01, double p11, double p21,
									  double D11, double D21, double p02,
									  double p12, double p22, double D12,
									  double D22, double[] xA, double[] xB,
									  int tau1, int tau2) {
		return quzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2);
	}

	private static double ppippispzderiv2(double p01, double p11, double p21,
										  double D11, double D21, double p02,
										  double p12, double p22, double D12,
										  double D22, double[] xA, double[] xB,
										  int tau1, int tau2) {
		return quzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) -
				uzQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
						xA, xB, tau1, tau2);
	}

	private static double pzpzspzderiv2(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double[] xA, double[] xB,
										int tau1, int tau2) {
		return quzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2) -
				uzQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
						xA, xB, tau1, tau2);
	}

	private static double sppisppideriv2(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return upiupideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
				xA, xB, tau1, tau2);
	}

	private static double spzspzderiv2(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double[] xA, double[] xB,
									   int tau1, int tau2) {
		return uzuzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA,
				xB, tau1, tau2);
	}

	private static double sppippipzderiv2(double p01, double p11, double p21,
										  double D11, double D21, double p02,
										  double p12, double p22, double D12,
										  double D22, double[] xA, double[] xB,
										  int tau1, int tau2) {
		return upiQpizderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
				xA, xB, tau1, tau2);
	}

	private static double ppipzsppideriv2(double p01, double p11, double p21,
										  double D11, double D21, double p02,
										  double p12, double p22, double D12,
										  double D22, double[] xA, double[] xB,
										  int tau1, int tau2) {
		return -upiQpizderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21,
				xA, xB, tau1, tau2);
	}

	private static double ppipzppipzderiv2(double p01, double p11, double p21,
										   double D11, double D21, double p02,
										   double p12, double p22, double D12,
										   double D22, double[] xA,
										   double[] xB,
										   int tau1, int tau2) {
		return QpizQpizderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
				xA, xB, tau1, tau2);
	}

	private static double pxpypxpyderiv2(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double[] xA, double[] xB,
										 int tau1, int tau2) {
		return 0.5 *
				(ppippippippideriv2(p01, p11, p21, D11, D21, p02, p12, p22,
						D12,
						D22, xA, xB, tau1, tau2) -
						pxpxpypyderiv2(p01, p11, p21, D11, D21, p02, p12, p22,
								D12, D22, xA, xB, tau1, tau2));
	}

	private static double LocalTwoCenterERIderiv2(NDDO6G a, NDDO6G b,
												  NDDO6G c, NDDO6G d,
												  int tau1, int tau2) {

		double[] A = a.getCoords();
		double[] C = c.getCoords();
		//(??|??)
		switch (a.getL()) {

			case 0://(s?|??)

				switch (b.getL()) {

					case 0: //(ss|??)

						switch (c.getL()) {

							case 0: //(ss|s?);

								switch (d.getL()) {

									case 0://(ss|ss)
										return ssssderiv2(a.p0, a.p1, a.p2,
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2, A, C, tau1, tau2);

									case 1:
										if (d.getk() == 1) {//(ss|spz)
											return ssspzderiv2(a.p0, a.p1,
													a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, A, C,
													tau1, tau2);
										}
										else {//(ss|sppi) = 0
											return 0;
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							case 1: //(ss|p?)
								if (c.getk() == 1) {//(ss|pz?)

									switch (d.getL()) {

										case 0://(ss|pzs)
											return ssspzderiv2(a.p0, a.p1,
													a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, A, C,
													tau1, tau2);

										case 1:
											if (d.getk() == 1) {//(ss|pzpz)
												return sspzpzderiv2(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														A, C, tau1, tau2);
											}
											else {//(ss|pzppi) = 0
												return 0;
											}
										default:
											return 0;
									}
								}
								else {//(ss|ppi?)

									if (d.getL() == 1 && d.getk() == 0 &&
											c.geti() == d.geti() &&
											c.getj() == d.getj()) {//(ss
										// |ppippi)
										return ssppippideriv2(a.p0, a.p1, a.p2,
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2, A, C, tau1, tau2);
									}
									else {//all others are 0
										return 0;
									}
								}
							default:
								System.err.println("oh no");
								return 0;

						}
					case 1: //(sp|??)

						if (b.getk() == 1) {//(spz|??)

							switch (c.getL()) {

								case 0://(spz|s?)

									switch (d.getL()) {

										case 0://(spz|ss)
											return spzssderiv2(a.p0, a.p1,
													a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, A, C,
													tau1, tau2);

										case 1:
											if (d.getk() == 1) {//(spz|spz)
												return spzspzderiv2(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														A, C, tau1, tau2);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(spz|pz?)

										switch (d.getL()) {

											case 0://(spz|pzs)
												return spzspzderiv2(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														A, C, tau1, tau2);

											case 1:
												if (d.getk() == 1) {//(spz
													// |pzpz)
													return spzpzpzderiv2(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, A, C, tau1,
															tau2);
												}
												else {//(spz|pzppi) = 0
													return 0;
												}
										}
									}
									else {//(spz|ppi?)
										if (d.geti() == c.geti() &&
												d.getj() == c.getj() &&
												d.getk() == 0) {
											return spzppippideriv2(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2, A,
													C, tau1, tau2);
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						}
						else {//(sppi|??)

							switch (c.getL()) {
								case 0://(sppi|s?)
									if (d.geti() == b.geti() &&
											d.getj() == b.getj() &&
											d.getk() == 0) {//(sppi|sppi)
										return sppisppideriv2(a.p0, a.p1, a.p2,
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2, A, C, tau1, tau2);
									}
									else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == b.geti() &&
												d.getj() == b.getj() &&
												d.getk() == 0) {//(sppi|pzppi)
											return sppippipzderiv2(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2, A,
													C, tau1, tau2);
										}
										else {
											return 0;
										}
									}
									else {
										if (c.geti() == b.geti() &&
												c.getj() == b.getj() &&
												c.getk() == 0) {//(sppi|ppi?)
											switch (d.getL()) {
												case 0:
													return sppisppideriv2(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, A, C, tau1,
															tau2);
												case 1:
													if (d.getk() == 1) {
														return sppippipzderiv2(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, A,
																C, tau1, tau2);
													}
													else {
														return 0;
													}
												default:
													return 0;
											}
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						}
					default:
						System.err.println("oh no");
						return 0;
				}

			case 1://(p?|??)
				if (a.getk() == 1) {//(pz?|??)
					switch (b.getL()) {
						case 0:
							switch (c.getL()) {

								case 0://(pzs|s?)

									switch (d.getL()) {

										case 0://(pzs|ss)
											return spzssderiv2(a.p0, a.p1,
													a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, A, C,
													tau1, tau2);

										case 1:
											if (d.getk() == 1) {//(pzs|spz)
												return spzspzderiv2(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														A, C, tau1, tau2);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(pzs|pz?)

										switch (d.getL()) {

											case 0://(pzs|pzs)
												return spzspzderiv2(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														A, C, tau1, tau2);

											case 1:
												if (d.getk() == 1) {//(pzs
													// |pzpz)
													return spzpzpzderiv2(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, A, C, tau1,
															tau2);
												}
												else {//(pzs|pzppi) = 0
													return 0;
												}
										}
									}
									else {//(pzs|ppi?)
										if (d.geti() == c.geti() &&
												d.getj() == c.getj() &&
												d.getk() == 0) {
											return spzppippideriv2(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2, A,
													C, tau1, tau2);
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:

							if (b.getk() == 1) {//(pzpz|??)

								switch (c.getL()) {

									case 0://(pzpz|s?)

										switch (d.getL()) {

											case 0://(pzpz|ss)
												return pzpzssderiv2(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														A, C, tau1, tau2);

											case 1:
												if (d.getk() == 1) {//(pzpz
													// |spz)
													return pzpzspzderiv2(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, A, C, tau1,
															tau2);
												}
												else {
													return 0;
												}
										}

									case 1:
										if (c.getk() == 1) {//(pzpz|pz?)

											switch (d.getL()) {

												case 0://(pzpz|pzs)
													return pzpzspzderiv2(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, A, C, tau1,
															tau2);

												case 1:
													if (d.getk() ==
															1) {//(pzpz|pzpz)
														return pzpzpzpzderiv2(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, A,
																C, tau1, tau2);
													}
													else {//(pzpz|pzppi) = 0
														return 0;
													}
											}
										}
										else {//(pzpz|ppi?)
											if (d.geti() == c.geti() &&
													d.getj() == c.getj() &&
													d.getk() == 0) {
												return pzpzppippideriv2(a.p0,
														a.p1, a.p2, a.D1, a.D2,
														c.p0, c.p1, c.p2, c.D1,
														c.D2, A, C, tau1,
														tau2);
											}
											else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							}
							else {//(pzppi|??)

								switch (c.getL()) {
									case 0://(pzppi|s?)
										if (d.geti() == b.geti() &&
												d.getj() == b.getj() &&
												d.getk() == 0) {//(pzppi|sppi)
											return ppipzsppideriv2(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2, A,
													C, tau1, tau2);
										}
										else {
											return 0;
										}
									case 1:
										if (c.getk() == 1) {
											if (d.geti() == b.geti() &&
													d.getj() == b.getj() &&
													d.getk() ==
															0) {//(pzppi|pzppi)
												return ppipzppipzderiv2(a.p0,
														a.p1, a.p2, a.D1, a.D2,
														c.p0, c.p1, c.p2, c.D1,
														c.D2, A, C, tau1,
														tau2);
											}
											else {
												return 0;
											}
										}
										else {
											if (c.geti() == b.geti() &&
													c.getj() == b.getj() &&
													c.getk() ==
															0) {//(pzppi|ppi?)
												switch (d.getL()) {
													case 0:
														return ppipzsppideriv2(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, A,
																C, tau1, tau2);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzderiv2(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	A, C, tau1,
																	tau2);
														}
														else {
															return 0;
														}
													default:
														return 0;
												}
											}
											else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							}
					}
				}
				else {//(ppi?|??);

					switch (b.getL()) {
						case 0://(ppis|??)

							switch (c.getL()) {
								case 0://(ppis|s?)
									if (d.geti() == a.geti() &&
											d.getj() == a.getj() &&
											d.getk() == 0) {//(ppis|sppi)
										return sppisppideriv2(a.p0, a.p1, a.p2,
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2, A, C, tau1, tau2);
									}
									else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == a.geti() &&
												d.getj() == a.getj() &&
												d.getk() == 0) {//(ppis|pzppi)
											return sppippipzderiv2(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2, A,
													C, tau1, tau2);
										}
										else {
											return 0;
										}
									}
									else {
										if (c.geti() == a.geti() &&
												c.getj() == a.getj() &&
												c.getk() == 0) {//(ppis|ppi?)
											switch (d.getL()) {
												case 0:
													return sppisppideriv2(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, A, C, tau1,
															tau2);
												case 1:
													if (d.getk() == 1) {
														return sppippipzderiv2(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, A,
																C, tau1, tau2);
													}
													else {
														return 0;
													}
												default:
													return 0;
											}
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:
							if (b.getk() == 1) {//(ppipz|??)
								switch (c.getL()) {
									case 0://(ppipz|s?)
										if (d.geti() == a.geti() &&
												d.getj() == a.getj() &&
												d.getk() == 0) {//(ppipz|sppi)
											return ppipzsppideriv2(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2, A,
													C, tau1, tau2);
										}
										else {
											return 0;
										}
									case 1:
										if (c.getk() == 1) {
											if (d.geti() == a.geti() &&
													d.getj() == a.getj() &&
													d.getk() ==
															0) {//(ppipz|pzppi)
												return ppipzppipzderiv2(a.p0,
														a.p1, a.p2, a.D1, a.D2,
														c.p0, c.p1, c.p2, c.D1,
														c.D2, A, C, tau1,
														tau2);
											}
											else {
												return 0;
											}
										}
										else {
											if (c.geti() == a.geti() &&
													c.getj() == a.getj() &&
													c.getk() ==
															0) {//(ppipz|ppi?)
												switch (d.getL()) {
													case 0:
														return ppipzsppideriv2(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, A,
																C, tau1, tau2);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzderiv2(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	A, C, tau1,
																	tau2);
														}
														else {
															return 0;
														}
													default:
														return 0;
												}
											}
											else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							}
							else {

								switch (c.getL()) {
									case 0://(ppippi|s?)
										switch (d.getL()) {
											case 0://(ppippi|ss)
												if (a.geti() == b.geti() &&
														a.getj() == b.getj()) {
													return ppippissderiv2(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, A, C, tau1,
															tau2);
												}
												else {
													return 0;
												}
											case 1:
												if (d.getk() == 1 &&
														a.geti() == b.geti() &&
														a.getj() == b.getj()) {
													return ppippispzderiv2(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, A, C, tau1,
															tau2);
												}
												else {
													return 0;
												}
										}

									case 1:
										if (c.getk() == 1) {
											switch (d.getL()) {
												case 0://(ppippi|pzs)
													if (a.geti() == b.geti() &&
															a.getj() ==
																	b.getj()) {
														return ppippispzderiv2(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, A,
																C, tau1, tau2);
													}
													else {
														return 0;
													}

												case 1:
													if (d.getk() == 1 &&
															a.geti() ==
																	b.geti() &&
															a.getj() ==
																	b.getj()) {//(ppippi|pzpz)
														return ppippipzpzderiv2(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, A,
																C, tau1, tau2);
													}
													else {
														return 0;
													}
											}
										}
										else {
											if (a.geti() == b.geti() &&
													a.getj() ==
															b.getj()) {//(pxpx
												// |??) or (pypy|??)

												if (c.getL() == d.getL() &&
														c.geti() == d.geti() &&
														c.getj() == d.getj() &&
														c.getk() == 0) {
													if (a.geti() == c.geti()) {
														return ppippippippideriv2(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, A,
																C, tau1, tau2);
													}
													else {
														return pxpxpypyderiv2(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, A,
																C, tau1, tau2);
													}
												}
												else {
													return 0;
												}

											}
											else {//(pxpy|??) or (pypx|??)
												if (c.getL() == d.getL() &&
														c.geti() != d.geti() &&
														c.getj() != d.getj() &&
														c.getk() == 0) {
													return pxpypxpyderiv2(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, A, C, tau1,
															tau2);
												}
											}
										}

								}

							}
					}

				}
		}

		return 0;
	}

	@Deprecated
	public static double getGderiv2finite(NDDO6G a, NDDO6G b, NDDO6G c,
										  NDDO6G d, int tau1, int tau2) {
		double orig = GeometryDerivative.getGderiv(a, b, c, d, tau1);

		double[] newcoords = a.getCoords().clone();

		newcoords[tau2] += 1E-6;

		NDDO6G anew = new NDDO6G(a, newcoords);
		NDDO6G bnew = new NDDO6G(b, newcoords);

		double perturbed = GeometryDerivative.getGderiv(anew, bnew, c, d,
				tau1);

		return (perturbed - orig) / 1E-6;
	}

	@Deprecated
	public static double[] secondDerivativeDecompositionfinite(double[] point1,
															   double[] point2,
															   NDDO6G a,
															   int tau1,
															   int tau2) {
		if (a.getL() == 0) {
			return new double[]{0};
		}

		double[] orig = GeometryDerivative
				.derivativeDecomposition(point1, point2, a, tau1);

		point1 = point1.clone();

		point1[tau2] += 1E-9;

		double[] perturbed = GeometryDerivative
				.derivativeDecomposition(point1, point2, a, tau1);

		return new double[]{(perturbed[0] - orig[0]) / 1E-9,
				(perturbed[1] - orig[1]) / 1E-9,
				(perturbed[2] - orig[2]) / 1E-9};
	}

	private static double[] secondDerivativeDecomposition(double[] point1,
														  double[] point2,
														  NDDO6G a, int tau1,
														  int tau2) {
		if (a.getL() == 0) {
			return new double[]{0};
		}

		int A = Math.min(tau1, tau2);
		int B = Math.max(tau1, tau2);

		double x = point2[0] - point1[0];
		double y = point2[1] - point1[1];
		double z = point2[2] - point1[2];

		double R = GTO.R(point1, point2);
		double Rxy = Math.sqrt(x * x + y * y);

		double[] returnval = new double[3];


		if (a.geti() == 1 && a.getL() == 1) {

			switch (A) {
				case 0:

					switch (B) {
						case 0://partial wrt x and x
							returnval[0] = 3 * x * z / (R * R * R * Rxy) +
									3 * x * z / (R * Rxy * Rxy * Rxy)
									- 2 * x * x * x * z /
									(R * R * R * Rxy * Rxy * Rxy)
									- 3 * x * x * x * z /
									(R * Rxy * Rxy * Rxy * Rxy * Rxy)
									- 3 * x * x * x * z /
									(R * R * R * R * R * Rxy);

							returnval[1] = y / (Rxy * Rxy * Rxy) -
									3 * x * x * y /
											(Rxy * Rxy * Rxy * Rxy * Rxy);

							returnval[2] =
									3 * x * x * x / (R * R * R * R * R) -
											3 * x / (R * R * R);


							return returnval;

						case 1://partial wrt x and y
							returnval[0] = z * y / (R * R * R * Rxy) +
									z * y / (R * Rxy * Rxy * Rxy)
									- 2 * x * x * y * z /
									(R * R * R * Rxy * Rxy * Rxy)
									- 3 * x * x * y * z /
									(R * R * R * R * R * Rxy)
									- 3 * x * x * y * z /
									(Rxy * Rxy * Rxy * Rxy * Rxy * R);

							returnval[1] = x / (Rxy * Rxy * Rxy) -
									3 * x * y * y /
											(Rxy * Rxy * Rxy * Rxy * Rxy);

							returnval[2] =
									3 * x * x * y / (R * R * R * R * R) -
											y / (R * R * R);


							return returnval;

						case 2://partial wrt x and z
							returnval[0] =
									z * z / (R * R * R * Rxy) - 1 / (R * Rxy)
											+ x * x / (R * R * R * Rxy) +
											x * x / (Rxy * Rxy * Rxy * R)
											- 3 * x * x * z * z /
											(R * R * R * R * R * Rxy)
											- x * x * z * z /
											(R * R * R * Rxy * Rxy * Rxy);

							returnval[1] = 0;

							returnval[2] =
									3 * x * x * z / (R * R * R * R * R) -
											z / (R * R * R);


							return returnval;

					}
					break;
				case 1:

					switch (B) {

						case 1://partial wrt y and y
							returnval[0] = x * z / (R * R * R * Rxy) +
									x * z / (R * Rxy * Rxy * Rxy)
									- 2 * x * y * y * z /
									(R * R * R * Rxy * Rxy * Rxy)
									- 3 * x * y * y * z /
									(R * R * R * R * R * Rxy)
									- 3 * x * y * y * z /
									(Rxy * Rxy * Rxy * Rxy * Rxy * R);

							returnval[1] = 3 * y / (Rxy * Rxy * Rxy) -
									3 * y * y * y /
											(Rxy * Rxy * Rxy * Rxy * Rxy);

							returnval[2] =
									3 * x * y * y / (R * R * R * R * R) -
											x / (R * R * R);

							return returnval;

						case 2://partial wrt y and z
							returnval[0] = x * y / (R * R * R * Rxy) +
									x * y / (R * Rxy * Rxy * Rxy)
									- 3 * x * y * z * z /
									(R * R * R * R * R * Rxy)
									- x * y * z * z /
									(R * R * R * Rxy * Rxy * Rxy);

							returnval[1] = 0;

							returnval[2] = 3 * x * y * z / (R * R * R * R * R);


							return returnval;

					}
				case 2:
					returnval[0] = 3 * x * z / (R * R * R * Rxy) -
							3 * x * z * z * z / (R * R * R * R * R * Rxy);

					returnval[1] = 0;

					returnval[2] = 3 * x * z * z / (R * R * R * R * R) -
							x / (R * R * R);


					return returnval;
				default:
			}
		}
		else if (a.getj() == 1 && a.getL() == 1) {

			switch (A) {
				case 0:

					switch (B) {
						case 0: //partial wrt x and x
							returnval[0] = y * z / (R * R * R * Rxy) +
									y * z / (R * Rxy * Rxy * Rxy)
									- 2 * x * x * y * z /
									(R * R * R * Rxy * Rxy * Rxy)
									- 3 * x * x * y * z /
									(R * R * R * R * R * Rxy)
									- 3 * x * x * y * z /
									(Rxy * Rxy * Rxy * Rxy * Rxy * R);

							returnval[1] = -3 * x / (Rxy * Rxy * Rxy) +
									3 * x * x * x /
											(Rxy * Rxy * Rxy * Rxy * Rxy);

							returnval[2] =
									3 * x * x * y / (R * R * R * R * R) -
											y / (R * R * R);

							return returnval;

						case 1://partial wrt x and y
							returnval[0] = z * x / (R * R * R * Rxy) +
									z * x / (R * Rxy * Rxy * Rxy)
									- 2 * x * y * y * z /
									(R * R * R * Rxy * Rxy * Rxy)
									- 3 * x * y * y * z /
									(R * R * R * R * R * Rxy)
									- 3 * x * y * y * z /
									(Rxy * Rxy * Rxy * Rxy * Rxy * R);

							returnval[1] = -y / (Rxy * Rxy * Rxy) +
									3 * x * x * y /
											(Rxy * Rxy * Rxy * Rxy * Rxy);

							returnval[2] =
									3 * x * y * y / (R * R * R * R * R) -
											x / (R * R * R);


							return returnval;

						case 2://partial wrt x and z

							returnval[0] = x * y / (R * R * R * Rxy) +
									x * y / (R * Rxy * Rxy * Rxy)
									- 3 * x * y * z * z /
									(R * R * R * R * R * Rxy)
									- x * y * z * z /
									(R * R * R * Rxy * Rxy * Rxy);

							returnval[1] = 0;

							returnval[2] = 3 * x * y * z / (R * R * R * R * R);


							return returnval;

					}
					break;
				case 1:

					switch (B) {

						case 1: //partial wrt y and y
							returnval[0] = 3 * y * z / (R * R * R * Rxy) +
									3 * y * z / (R * Rxy * Rxy * Rxy)
									- 2 * y * y * y * z /
									(R * R * R * Rxy * Rxy * Rxy)
									- 3 * y * y * y * z /
									(R * Rxy * Rxy * Rxy * Rxy * Rxy)
									- 3 * y * y * y * z /
									(R * R * R * R * R * Rxy);

							returnval[1] = -x / (Rxy * Rxy * Rxy) +
									3 * x * y * y /
											(Rxy * Rxy * Rxy * Rxy * Rxy);

							returnval[2] =
									3 * y * y * y / (R * R * R * R * R) -
											3 * y / (R * R * R);

							return returnval;

						case 2://partial wrt y and z
							returnval[0] =
									z * z / (R * R * R * Rxy) - 1 / (R * Rxy)
											+ y * y / (R * R * R * Rxy) +
											y * y / (Rxy * Rxy * Rxy * R)
											- 3 * y * y * z * z /
											(R * R * R * R * R * Rxy)
											- y * y * z * z /
											(R * R * R * Rxy * Rxy * Rxy);

							returnval[1] = 0;

							returnval[2] =
									3 * y * y * z / (R * R * R * R * R) -
											z / (R * R * R);


							return returnval;

					}
				case 2:
					returnval[0] = 3 * y * z / (R * R * R * Rxy) -
							3 * y * z * z * z / (R * R * R * R * R * Rxy);

					returnval[1] = 0;

					returnval[2] = 3 * y * z * z / (R * R * R * R * R) -
							y / (R * R * R);


					return returnval;
				default:
			}
		}
		else if (a.getk() == 1 && a.getL() == 1) {

			switch (A) {
				case 0:

					switch (B) {

						case 0: //partial wrt x and x
							returnval[0] =
									3 * Rxy * x * x / (R * R * R * R * R) -
											Rxy / (R * R * R) + 1 / (Rxy * R)
											- 2 * x * x / (Rxy * R * R * R) -
											x * x / (R * Rxy * Rxy * Rxy);

							returnval[1] = 0;

							returnval[2] =
									3 * x * x * z / (R * R * R * R * R) -
											z / (R * R * R);

							return returnval;

						case 1://partial wrt x and y
							returnval[0] =
									3 * Rxy * x * y / (R * R * R * R * R) -
											2 * x * y / (Rxy * R * R * R)
											- x * y / (R * Rxy * Rxy * Rxy);

							returnval[1] = 0;

							returnval[2] = 3 * x * y * z / (R * R * R * R * R);


							return returnval;

						case 2://partial wrt x and z

							returnval[0] =
									3 * Rxy * x * z / (R * R * R * R * R) -
											x * z / (Rxy * R * R * R);

							returnval[1] = 0;

							returnval[2] =
									3 * x * z * z / (R * R * R * R * R) -
											x / (R * R * R);


							return returnval;

					}
					break;
				case 1:

					switch (B) {

						case 1: //partial wrt y and y
							returnval[0] =
									3 * Rxy * y * y / (R * R * R * R * R) -
											Rxy / (R * R * R) + 1 / (Rxy * R)
											- 2 * y * y / (Rxy * R * R * R) -
											y * y / (R * Rxy * Rxy * Rxy);

							returnval[1] = 0;

							returnval[2] =
									3 * y * y * z / (R * R * R * R * R) -
											z / (R * R * R);

							return returnval;

						case 2://partial wrt y and z

							returnval[0] =
									3 * Rxy * y * z / (R * R * R * R * R) -
											y * z / (Rxy * R * R * R);

							returnval[1] = 0;

							returnval[2] =
									3 * y * z * z / (R * R * R * R * R) -
											y / (R * R * R);


							return returnval;

					}
				case 2:
					returnval[0] = 3 * z * z * Rxy / (R * R * R * R * R) -
							Rxy / (R * R * R);

					returnval[1] = 0;

					returnval[2] = 3 * z * z * z / (R * R * R * R * R) -
							3 * z / (R * R * R);

					return returnval;
				default:
			}
		}

		System.err.println("oh no!");

		return new double[]{0, 0, 0};


	}

	@Deprecated
	public static double[] secondDerivativeDecomposition2finite(double[] point1,
																double[] point2,
																NDDO6G a,
																int tau1,
																int tau2) {
		if (a.getL() == 0) {
			return new double[]{0};
		}

		double[] orig = GeometryDerivative
				.derivativeDecomposition2(point1, point2, a, tau1);

		point1 = point1.clone();

		point1[tau2] += 1E-9;

		double[] perturbed = GeometryDerivative
				.derivativeDecomposition2(point1, point2, a, tau1);

		return new double[]{(perturbed[0] - orig[0]) / 1E-9,
				(perturbed[1] - orig[1]) / 1E-9,
				(perturbed[2] - orig[2]) / 1E-9};
	}

	private static double[] secondDerivativeDecomposition2(double[] point1,
														   double[] point2,
														   NDDO6G a, int tau1,
														   int tau2) {
		if (a.getL() == 0) {
			return new double[]{0};
		}

		int A = Math.min(tau1, tau2);
		int B = Math.max(tau1, tau2);

		double x = point2[0] - point1[0];
		double y = point2[1] - point1[1];
		double z = point2[2] - point1[2];

		double R = GTO.R(point1, point2);
		double Rxz = Math.sqrt(x * x + z * z);

		double[] returnval = new double[3];

		if (a.geti() == 1 && a.getL() == 1) {

			switch (A) {
				case 0:

					switch (B) {
						case 0://partial wrt x and x
							returnval[0] = 3 * x * x * x * y /
									(R * Rxz * Rxz * Rxz * Rxz * Rxz)
									+ 3 * x * x * x * y /
									(R * R * R * R * R * Rxz)
									+ 2 * x * x * x * y /
									(R * R * R * Rxz * Rxz * Rxz)
									- 3 * x * y / (R * R * R * Rxz)
									- 3 * x * y / (R * Rxz * Rxz * Rxz);

							returnval[1] = z / (Rxz * Rxz * Rxz) -
									3 * x * x * z /
											(Rxz * Rxz * Rxz * Rxz * Rxz);

							returnval[2] =
									3 * x * x * x / (R * R * R * R * R) -
											3 * x / (R * R * R);


							return returnval;

						case 1://partial wrt x and y
							returnval[0] = 1 / (R * Rxz) + 3 * x * x * y * y /
									(R * R * R * R * R * Rxz)
									+ x * x * y * y /
									(R * R * R * Rxz * Rxz * Rxz)
									- x * x / (R * Rxz * Rxz * Rxz) -
									x * x / (R * R * R * Rxz)
									- y * y / (R * R * R * Rxz);

							returnval[1] = 0;

							returnval[2] =
									3 * x * x * y / (R * R * R * R * R) -
											y / (R * R * R);


							return returnval;

						case 2://partial wrt x and z
							returnval[0] = 3 * x * x * y * z /
									(R * R * R * R * R * Rxz)
									+ 2 * x * x * y * z /
									(R * R * R * Rxz * Rxz * Rxz)
									+ 3 * x * x * y * z /
									(R * Rxz * Rxz * Rxz * Rxz * Rxz)
									- y * z / (R * R * R * Rxz)
									- y * z / (R * Rxz * Rxz * Rxz);


							returnval[1] = x / (Rxz * Rxz * Rxz) -
									3 * x * z * z /
											(Rxz * Rxz * Rxz * Rxz * Rxz);

							returnval[2] =
									3 * x * x * z / (R * R * R * R * R) -
											z / (R * R * R);


							return returnval;

					}
					break;
				case 1:

					switch (B) {

						case 1://partial wrt y and y
							returnval[0] = 3 * x * y * y * y /
									(R * R * R * R * R * Rxz)
									- 3 * x * y / (R * R * R * Rxz);

							returnval[1] = 0;

							returnval[2] =
									3 * x * y * y / (R * R * R * R * R) -
											x / (R * R * R);

							return returnval;

						case 2://partial wrt y and z
							returnval[0] = 3 * x * y * y * z /
									(R * R * R * R * R * Rxz)
									+ x * y * y * z /
									(R * R * R * Rxz * Rxz * Rxz)
									- x * z / (R * Rxz * Rxz * Rxz)
									- x * z / (R * R * R * Rxz);

							returnval[1] = 0;

							returnval[2] = 3 * x * y * z / (R * R * R * R * R);


							return returnval;

					}
				case 2:
					returnval[0] =
							3 * x * y * z * z / (R * R * R * R * R * Rxz)
									+ 3 * x * y * z * z /
									(R * Rxz * Rxz * Rxz * Rxz * Rxz)
									+ 2 * x * y * z * z /
									(R * R * R * Rxz * Rxz * Rxz)
									- x * y / (R * Rxz * Rxz * Rxz) -
									x * y / (R * R * R * Rxz);

					returnval[1] = 3 * z / (Rxz * Rxz * Rxz) -
							3 * z * z * z / (Rxz * Rxz * Rxz * Rxz * Rxz);

					returnval[2] = 3 * x * z * z / (R * R * R * R * R) -
							x / (R * R * R);


					return returnval;
				default:
			}
		}
		else if (a.getj() == 1 && a.getL() == 1) {

			switch (A) {
				case 0:

					switch (B) {
						case 0: //partial wrt x and x
							returnval[0] = 2 * x * x / (R * R * R * Rxz) +
									x * x / (R * Rxz * Rxz * Rxz)
									+ Rxz / (R * R * R) - 1 / (R * Rxz) -
									3 * x * x * Rxz / (R * R * R * R * R);
							returnval[1] = 0;

							returnval[2] =
									3 * x * x * y / (R * R * R * R * R) -
											y / (R * R * R);

							return returnval;

						case 1://partial wrt x and y
							returnval[0] = x * y / (R * R * R * Rxz) -
									3 * x * y * Rxz / (R * R * R * R * R);

							returnval[1] = 0;

							returnval[2] =
									3 * x * y * y / (R * R * R * R * R) -
											x / (R * R * R);


							return returnval;

						case 2://partial wrt x and z

							returnval[0] = 2 * x * z / (R * R * R * Rxz) +
									x * z / (R * Rxz * Rxz * Rxz)
									- 3 * x * z * Rxz / (R * R * R * R * R);

							returnval[1] = 0;

							returnval[2] = 3 * x * y * z / (R * R * R * R * R);


							return returnval;

					}
					break;
				case 1:

					switch (B) {

						case 1: //partial wrt y and y
							returnval[0] = Rxz / (R * R * R) -
									3 * y * y * Rxz / (R * R * R * R * R);

							returnval[1] = 0;

							returnval[2] =
									3 * y * y * y / (R * R * R * R * R) -
											3 * y / (R * R * R);

							return returnval;

						case 2://partial wrt y and z
							returnval[0] = y * z / (R * R * R * Rxz) -
									3 * y * z * Rxz / (R * R * R * R * R);

							returnval[1] = 0;

							returnval[2] =
									3 * y * y * z / (R * R * R * R * R) -
											z / (R * R * R);


							return returnval;

					}
				case 2:
					returnval[0] = 2 * z * z / (R * R * R * Rxz) +
							z * z / (R * Rxz * Rxz * Rxz)
							+ Rxz / (R * R * R) - 1 / (R * Rxz) -
							3 * z * z * Rxz / (R * R * R * R * R);

					returnval[1] = 0;

					returnval[2] = 3 * y * z * z / (R * R * R * R * R) -
							y / (R * R * R);


					return returnval;
				default:
			}
		}
		else if (a.getk() == 1 && a.getL() == 1) {

			switch (A) {
				case 0:

					switch (B) {

						case 0: //partial wrt x and x
							returnval[0] = 2 * x * x * y * z /
									(R * R * R * Rxz * Rxz * Rxz) +
									3 * x * x * y * z /
											(R * R * R * R * R * Rxz)
									+ 3 * x * x * y * z /
									(R * Rxz * Rxz * Rxz * Rxz * Rxz) -
									y * z / (R * R * R * Rxz) -
									y * z / (R * Rxz * Rxz * Rxz);

							returnval[1] = 3 * x * x * x /
									(Rxz * Rxz * Rxz * Rxz * Rxz) -
									3 * x / (Rxz * Rxz * Rxz);

							returnval[2] =
									3 * x * x * z / (R * R * R * R * R) -
											z / (R * R * R);

							return returnval;

						case 1://partial wrt x and y
							returnval[0] = 3 * x * y * y * z /
									(R * R * R * R * R * Rxz)
									+ x * y * y * z /
									(R * R * R * Rxz * Rxz * Rxz)
									- x * z / (R * Rxz * Rxz * Rxz) -
									x * z / (R * R * R * Rxz);

							returnval[1] = 0;

							returnval[2] = 3 * x * y * z / (R * R * R * R * R);


							return returnval;

						case 2://partial wrt x and z

							returnval[0] = 2 * x * y * z * z /
									(R * R * R * Rxz * Rxz * Rxz) +
									3 * x * y * z * z /
											(R * R * R * R * R * Rxz)
									+ 3 * x * y * z * z /
									(R * Rxz * Rxz * Rxz * Rxz * Rxz) -
									x * y / (R * Rxz * Rxz * Rxz) -
									x * y / (R * R * R * Rxz);

							returnval[1] = 3 * x * x * z /
									(Rxz * Rxz * Rxz * Rxz * Rxz) -
									z / (Rxz * Rxz * Rxz);

							returnval[2] =
									3 * x * z * z / (R * R * R * R * R) -
											x / (R * R * R);


							return returnval;

					}
					break;
				case 1:

					switch (B) {

						case 1: //partial wrt y and y
							returnval[0] = 3 * y * y * y * z /
									(R * R * R * R * R * Rxz) -
									3 * y * z / (R * R * R * Rxz);

							returnval[1] = 0;

							returnval[2] =
									3 * y * y * z / (R * R * R * R * R) -
											z / (R * R * R);

							return returnval;

						case 2://partial wrt y and z

							returnval[0] = 1 / (R * Rxz) + 3 * y * y * z * z /
									(R * R * R * R * R * Rxz)
									+ y * y * z * z /
									(R * R * R * Rxz * Rxz * Rxz) -
									y * y / (R * R * R * Rxz)
									- z * z / (R * Rxz * Rxz * Rxz) -
									z * z / (R * R * R * Rxz);

							returnval[1] = 0;

							returnval[2] =
									3 * y * z * z / (R * R * R * R * R) -
											y / (R * R * R);


							return returnval;

					}
				case 2:
					returnval[0] =
							3 * y * z * z * z / (R * R * R * R * R * Rxz) +
									2 * y * z * z * z /
											(R * R * R * Rxz * Rxz * Rxz)
									+ 3 * y * z * z * z /
									(R * Rxz * Rxz * Rxz * Rxz * Rxz) -
									3 * y * z / (R * Rxz * Rxz * Rxz) -
									3 * y * z / (R * R * R * Rxz);

					returnval[1] =
							3 * x * z * z / (Rxz * Rxz * Rxz * Rxz * Rxz) -
									x / (Rxz * Rxz * Rxz);

					returnval[2] = 3 * z * z * z / (R * R * R * R * R) -
							3 * z / (R * R * R);

					return returnval;
				default:
			}
		}

		System.err.println("oh no!");

		return new double[]{0, 0, 0};


	}

	public static double getGderiv2(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d,
									int tau1, int tau2) {


		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double[] coeffAderiv1 = GeometryDerivative
				.derivativeDecomposition(a.getCoords(), c.getCoords(), a,
						tau1);
		double[] coeffBderiv1 = GeometryDerivative
				.derivativeDecomposition(a.getCoords(), c.getCoords(), b,
						tau1);
		double[] coeffCderiv1 = GeometryDerivative
				.derivativeDecomposition(a.getCoords(), c.getCoords(), c,
						tau1);
		double[] coeffDderiv1 = GeometryDerivative
				.derivativeDecomposition(a.getCoords(), c.getCoords(), d,
						tau1);

		double[] coeffAderiv2 = GeometryDerivative
				.derivativeDecomposition(a.getCoords(), c.getCoords(), a,
						tau2);
		double[] coeffBderiv2 = GeometryDerivative
				.derivativeDecomposition(a.getCoords(), c.getCoords(), b,
						tau2);
		double[] coeffCderiv2 = GeometryDerivative
				.derivativeDecomposition(a.getCoords(), c.getCoords(), c,
						tau2);
		double[] coeffDderiv2 = GeometryDerivative
				.derivativeDecomposition(a.getCoords(), c.getCoords(), d,
						tau2);

		double[] coeffAderiv =
				secondDerivativeDecomposition(a.getCoords(), c.getCoords(), a,
						tau1, tau2);
		double[] coeffBderiv =
				secondDerivativeDecomposition(a.getCoords(), c.getCoords(), b,
						tau1, tau2);
		double[] coeffCderiv =
				secondDerivativeDecomposition(a.getCoords(), c.getCoords(), c,
						tau1, tau2);
		double[] coeffDderiv =
				secondDerivativeDecomposition(a.getCoords(), c.getCoords(), d,
						tau1, tau2);

		if (Math.abs(a.getCoords()[0] - c.getCoords()[0]) < 1E-3 &&
				Math.abs(a.getCoords()[1] - c.getCoords()[1]) < 1E-3) {

			coeffAderiv =
					secondDerivativeDecomposition2(a.getCoords(),
							c.getCoords(),
							a, tau1, tau2);
			coeffBderiv =
					secondDerivativeDecomposition2(a.getCoords(),
							c.getCoords(),
							b, tau1, tau2);
			coeffCderiv =
					secondDerivativeDecomposition2(a.getCoords(),
							c.getCoords(),
							c, tau1, tau2);
			coeffDderiv =
					secondDerivativeDecomposition2(a.getCoords(),
							c.getCoords(),
							d, tau1, tau2);

			coeffAderiv2 = GeometryDerivative
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), a,
							tau2);
			coeffBderiv2 = GeometryDerivative
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), b,
							tau2);
			coeffCderiv2 = GeometryDerivative
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), c,
							tau2);
			coeffDderiv2 = GeometryDerivative
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), d,
							tau2);

			coeffAderiv1 = GeometryDerivative
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), a,
							tau1);
			coeffBderiv1 = GeometryDerivative
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), b,
							tau1);
			coeffCderiv1 = GeometryDerivative
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), c,
							tau1);
			coeffDderiv1 = GeometryDerivative
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), d,
							tau1);

			coeffA = a.decomposition2(a.getCoords(), c.getCoords());
			coeffB = b.decomposition2(a.getCoords(), c.getCoords());
			coeffC = c.decomposition2(a.getCoords(), c.getCoords());
			coeffD = d.decomposition2(a.getCoords(), c.getCoords());


		}


		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		double sum = 0;

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {

						double eri = NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k],
								D[l]);

						double erideriv1 = GeometryDerivative
								.LocalTwoCenterERIderiv(A[i], B[j], C[k], D[l],
										tau1);

						double erideriv2 = GeometryDerivative
								.LocalTwoCenterERIderiv(A[i], B[j], C[k], D[l],
										tau2);

						double erideriv =
								LocalTwoCenterERIderiv2(A[i], B[j], C[k], D[l],
										tau1, tau2);

						sum += coeffAderiv[i] * coeffB[j] * coeffC[k] *
								coeffD[l] * eri * 27.21;
						sum += coeffAderiv1[i] * coeffBderiv2[j] * coeffC[k] *
								coeffD[l] * eri * 27.21;
						sum += coeffAderiv1[i] * coeffB[j] * coeffCderiv2[k] *
								coeffD[l] * eri * 27.21;
						sum += coeffAderiv1[i] * coeffB[j] * coeffC[k] *
								coeffDderiv2[l] * eri * 27.21;
						sum += coeffAderiv1[i] * coeffB[j] * coeffC[k] *
								coeffD[l] * erideriv2 * 27.21;

						sum += coeffAderiv2[i] * coeffBderiv1[j] * coeffC[k] *
								coeffD[l] * eri * 27.21;
						sum += coeffA[i] * coeffBderiv[j] * coeffC[k] *
								coeffD[l] * eri * 27.21;
						sum += coeffA[i] * coeffBderiv1[j] * coeffCderiv2[k] *
								coeffD[l] * eri * 27.21;
						sum += coeffA[i] * coeffBderiv1[j] * coeffC[k] *
								coeffDderiv2[l] * eri * 27.21;
						sum += coeffA[i] * coeffBderiv1[j] * coeffC[k] *
								coeffD[l] * erideriv2 * 27.21;

						sum += coeffAderiv2[i] * coeffB[j] * coeffCderiv1[k] *
								coeffD[l] * eri * 27.21;
						sum += coeffA[i] * coeffBderiv2[j] * coeffCderiv1[k] *
								coeffD[l] * eri * 27.21;
						sum += coeffA[i] * coeffB[j] * coeffCderiv[k] *
								coeffD[l] * eri * 27.21;
						sum += coeffA[i] * coeffB[j] * coeffCderiv1[k] *
								coeffDderiv2[l] * eri * 27.21;
						sum += coeffA[i] * coeffB[j] * coeffCderiv1[k] *
								coeffD[l] * erideriv2 * 27.21;

						sum += coeffAderiv2[i] * coeffB[j] * coeffC[k] *
								coeffDderiv1[l] * eri * 27.21;
						sum += coeffA[i] * coeffBderiv2[j] * coeffC[k] *
								coeffDderiv1[l] * eri * 27.21;
						sum += coeffA[i] * coeffB[j] * coeffCderiv2[k] *
								coeffDderiv1[l] * eri * 27.21;
						sum += coeffA[i] * coeffB[j] * coeffC[k] *
								coeffDderiv[l] * eri * 27.21;
						sum += coeffA[i] * coeffB[j] * coeffC[k] *
								coeffDderiv1[l] * erideriv2 * 27.21;

						sum += coeffAderiv2[i] * coeffB[j] * coeffC[k] *
								coeffD[l] * erideriv1 * 27.21;
						sum += coeffA[i] * coeffBderiv2[j] * coeffC[k] *
								coeffD[l] * erideriv1 * 27.21;
						sum += coeffA[i] * coeffB[j] * coeffCderiv2[k] *
								coeffD[l] * erideriv1 * 27.21;
						sum += coeffA[i] * coeffB[j] * coeffC[k] *
								coeffDderiv2[l] * erideriv1 * 27.21;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
								erideriv * 27.21;


					}
				}
			}
		}

		return sum;
	}

	private static double Ederiv2(int atomnum1, int atomnum2, int[][] index,
								  SimpleMatrix densityMatrix, NDDOAtom[] atoms,
								  NDDO6G[] orbitals, int tau1, int tau2) {

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += densityMatrix.get(i, j) * atoms[atomnum2]
							.Vderiv2(orbitals[i], orbitals[j], tau1, tau2);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e += densityMatrix.get(k, l) * atoms[atomnum1]
							.Vderiv2(orbitals[k], orbitals[l], tau1, tau2);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * densityMatrix.get(i, k) *
							NDDO6G.betaderiv2(orbitals[i], orbitals[k], tau1,
									tau2);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {
							e += (densityMatrix.get(i, j) *
									densityMatrix.get(k, l) -
									densityMatrix.get(i, k) * 0.5 *
											densityMatrix.get(j, l))
									* GeometrySecondDerivative
									.getGderiv2(orbitals[i], orbitals[j],
											orbitals[k], orbitals[l], tau1,
											tau2);
						}
					}
				}
			}
		}

		return e;


	}

	private static double Ederiv2(int atomnum1, int atomnum2, int[][] index,
								  SimpleMatrix alphaDensity,
								  SimpleMatrix betaDensity, NDDOAtom[] atoms,
								  NDDO6G[] orbitals, int tau1, int tau2) {

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += (alphaDensity.get(i, j) + betaDensity.get(i, j)) *
							atoms[atomnum2]
									.Vderiv2(orbitals[i], orbitals[j], tau1,
											tau2);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e += (alphaDensity.get(k, l) + betaDensity.get(k, l)) *
							atoms[atomnum1]
									.Vderiv2(orbitals[k], orbitals[l], tau1,
											tau2);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * (alphaDensity.get(i, k) + betaDensity.get(i, k)) *
							NDDO6G.betaderiv2(orbitals[i], orbitals[k], tau1,
									tau2);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {
							e += ((alphaDensity.get(i, j) +
									betaDensity.get(i, j)) *
									(alphaDensity.get(k, l) +
											betaDensity.get(k, l)) -
									alphaDensity.get(i, k) *
											alphaDensity.get(j, l) -
									betaDensity.get(i, k) *
											betaDensity.get(j, l))
									* GeometrySecondDerivative
									.getGderiv2(orbitals[i], orbitals[j],
											orbitals[k], orbitals[l], tau1,
											tau2);
						}
					}
				}
			}
		}

		return e;


	}

	@Deprecated
	public static double hessianfinite(NDDOAtom[] atoms, SolutionR soln,
									   int atomnum1, int tau1, int atomnum2,
									   int tau2) {

		double initval =
				GeometryDerivative.gradient(atoms, soln, atomnum1, tau1);

		NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum2, tau2);

		SolutionR newsoln = (SolutionR) soln.withNewAtoms(newatoms);

		double finalval =
				GeometryDerivative.gradient(newatoms, newsoln, atomnum1, tau1);

		return 1E7 * (finalval - initval);


	}

	@Deprecated
	public static double hessianfinite(NDDOAtom[] atoms, SolutionU soln,
									   int atomnum1, int tau1, int atomnum2,
									   int tau2) {

		double initval = GeometryDerivative
				.gradientUnrestricted(atoms, soln, atomnum1, tau1);

		NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum2, tau2);

		SolutionU newsoln = (SolutionU) soln.withNewAtoms(newatoms);

		double finalval = GeometryDerivative
				.gradientUnrestricted(newatoms, newsoln, atomnum1, tau1);

		return 1E7 * (finalval - initval);


	}

	public static SimpleMatrix hessianRoutine(SolutionR soln,
											  SimpleMatrix[] fockDerivStatic) {
		@SuppressWarnings("DuplicatedCode")
		SimpleMatrix[] densityDerivs =
				new SimpleMatrix[fockDerivStatic.length];
		int elapsedSize = 0;
		double cores = Runtime.getRuntime().availableProcessors();
		int size = Math.max((int) Math.ceil(fockDerivStatic.length / cores),
				1);

		List<RecursiveAction> subtasks = new ArrayList<>();
		// partitions densityDerivs into batches.
		while (elapsedSize < fockDerivStatic.length) {
			int finalElapsedSize = elapsedSize;
			subtasks.add(new RecursiveAction() {
				@Override
				protected void compute() {
					SimpleMatrix[] subset = Arrays.copyOfRange(fockDerivStatic,
							finalElapsedSize, Math.min(fockDerivStatic.length,
									finalElapsedSize + size));

					SimpleMatrix[] output;
					try {
						output = densityDerivPople(soln, subset);
					} catch (SingularMatrixException e) {
						output = densityDerivThiel(soln, subset);
					}

					// removed .dup() here
					System.arraycopy(output, 0, densityDerivs,
							finalElapsedSize, output.length);
				}
			});
			elapsedSize += size;
		}
		ForkJoinTask.invokeAll(subtasks);

		SimpleMatrix hessian =
				new SimpleMatrix(densityDerivs.length, densityDerivs.length);
		for (int i = 0; i < hessian.numRows(); i++) {
			for (int j = i; j < hessian.numRows(); j++) {
				double E = 0;
				int atomnum1 = i / 3;
				int atomnum2 = j / 3;
				int tau1 = i - 3 * atomnum1;
				int tau2 = j - 3 * atomnum2;

				if (atomnum1 == atomnum2) {
					for (int a = 0; a < soln.atoms.length; a++) {
						if (a != atomnum1) {
							E += Ederiv2(atomnum1, a, soln.orbsOfAtom,
									soln.densityMatrix(), soln.atoms,
									soln.orbitals,
									tau1, tau2);
							E += soln.atoms[atomnum1]
									.crfDeriv2(soln.atoms[a], tau1, tau2);
						}
					}
				}
				else {
					E = -Ederiv2(atomnum1, atomnum2, soln.orbsOfAtom,
							soln.densityMatrix(), soln.atoms, soln.orbitals,
							tau1,
							tau2) - soln.atoms[atomnum1]
							.crfDeriv2(soln.atoms[atomnum2], tau1, tau2);

				}

				for (int I = 0; I < soln.orbitals.length; I++) {
					for (int J = 0; J < soln.orbitals.length; J++) {
						E += fockDerivStatic[i].get(I, J) *
								densityDerivs[j].get(I, J);
					}
				}

				hessian.set(i, j, E);
				hessian.set(j, i, E);
			}
		}
		return hessian;
	}

	public static SimpleMatrix hessianRoutine(SolutionU soln,
											  SimpleMatrix[] fockderivstaticalpha,
											  SimpleMatrix[] fockderivstaticbeta) {


		SimpleMatrix[] densityderivsalpha =
				new SimpleMatrix[fockderivstaticalpha.length];
		SimpleMatrix[] densityderivsbeta =
				new SimpleMatrix[fockderivstaticbeta.length];


		int count = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			for (int tau = 0; tau < 3; tau++) {

				SimpleMatrix[] matrices = GeometryDerivative
						.densitymatrixderivfinite(soln.atoms, soln, a, tau);
				densityderivsalpha[count] = matrices[0];
				densityderivsbeta[count] = matrices[1];
				count++;
			}
		}

		SimpleMatrix hessian = new SimpleMatrix(densityderivsalpha.length,
				densityderivsalpha.length);

		for (int i = 0; i < hessian.numRows(); i++) {
			for (int j = i; j < hessian.numRows(); j++) {

				double E = 0;

				int atomnum1 = i / 3;

				int atomnum2 = j / 3;

				int tau1 = i - 3 * atomnum1;

				int tau2 = j - 3 * atomnum2;

				if (atomnum1 == atomnum2) {
					for (int a = 0; a < soln.atoms.length; a++) {
						if (a != atomnum1) {
							E += Ederiv2(atomnum1, a, soln.orbsOfAtom,
									soln.alphaDensity(), soln.betaDensity(),
									soln.atoms, soln.orbitals, tau1, tau2);
							E += soln.atoms[atomnum1]
									.crfDeriv2(soln.atoms[a], tau1, tau2);
						}
					}
				}
				else {
					E = -Ederiv2(atomnum1, atomnum2, soln.orbsOfAtom,
							soln.alphaDensity(), soln.betaDensity(),
							soln.atoms,
							soln.orbitals, tau1, tau2) - soln.atoms[atomnum1]
							.crfDeriv2(soln.atoms[atomnum2], tau1, tau2);

				}


				for (int I = 0; I < soln.orbitals.length; I++) {
					for (int J = 0; J < soln.orbitals.length; J++) {
						E += fockderivstaticalpha[i].get(I, J) *
								densityderivsalpha[j].get(I, J);
						E += fockderivstaticbeta[i].get(I, J) *
								densityderivsbeta[j].get(I, J);
					}
				}

				hessian.set(i, j, E);
				hessian.set(j, i, E);
			}
		}

		return hessian;


	}

	public static SimpleMatrix[] densityDerivThiel(SolutionR soln,
												   SimpleMatrix[] fockderivstatic) {
		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix[] xarray = new SimpleMatrix[fockderivstatic.length];
		SimpleMatrix[] rarray = new SimpleMatrix[fockderivstatic.length];
		SimpleMatrix[] dirs = new SimpleMatrix[fockderivstatic.length];

		double[] arrpreconditioner = new double[NOcc * NVirt];

		int counter = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = -soln.E.get(i) - soln.E.get(NOcc + j);

				arrpreconditioner[counter] = Math.pow(e, -0.5);

				counter++;
			}
		}

		SimpleMatrix D = SimpleMatrix.diag(arrpreconditioner);

		for (int a = 0; a < xarray.length; a++) {
			SimpleMatrix F = new SimpleMatrix(NOcc * NVirt, 1);

			int count1 = 0;

			for (int i = 0; i < NOcc; i++) {
				for (int j = 0; j < NVirt; j++) {

					double element = 0;

					for (int u = 0; u < soln.orbitals.length; u++) {
						for (int v = 0; v < soln.orbitals.length; v++) {
							element += soln.C.get(i, u) * soln.C.get(j + NOcc,
									v) * fockderivstatic[a].get(u, v);
						}
					}
					F.set(count1, 0, element);

					count1++;
				}
			}

			F = D.mult(F);

			xarray[a] = new SimpleMatrix(NOcc * NVirt, 1);
			rarray[a] = F;
			dirs[a] = F;
		}


		if (dirs[0].numRows() == 0) {
			SimpleMatrix[] densityderivs =
					new SimpleMatrix[fockderivstatic.length];

			for (int i = 0; i < densityderivs.length; i++) {
				densityderivs[i] = new SimpleMatrix(0, 0);
			}

			return densityderivs;
		}


		while (Utils.numNotNull(rarray) > 0) {
			ArrayList<SimpleMatrix> d = new ArrayList<>();

			ArrayList<SimpleMatrix> p = new ArrayList<>();

			for (int i = 0; i < rarray.length; i++) {
				if (rarray[i] != null) {
					d.add(new SimpleMatrix(dirs[i]));
					p.add(D.mult(
							computeResponseVectorsThiel(dirs[i], soln)));
				}
			}

			SimpleMatrix solver =
					new SimpleMatrix(p.size(), p.size());
			SimpleMatrix rhsvec =
					new SimpleMatrix(p.size(), rarray.length);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					double[] arrrhs = new double[p.size()];

					for (int i = 0; i < arrrhs.length; i++) {
						arrrhs[i] = 2 *
								rarray[a].transpose().mult(d.get(i))
										.get(0, 0);

					}
					rhsvec.setColumn(a, 0, arrrhs);
				}
			}

			for (int i = 0; i < solver.numRows(); i++) {
				for (int j = i; j < solver.numRows(); j++) {
					double val2 =
							p.get(j).transpose().mult(d.get(i))
									.get(0, 0) + p.get(i).transpose()
									.mult(d.get(j)).get(0, 0);

					solver.set(i, j, val2);
					solver.set(j, i, val2);
				}
			}

			SimpleMatrix alpha;
			try {
				alpha = solver.solve(rhsvec);
			} catch (SingularMatrixException e) {
				alpha = SimpleMatrix.ones(solver.numCols(), rhsvec.numCols());
			}

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					for (int i = 0; i < alpha.numRows(); i++) {
						xarray[a] = xarray[a].plus(
								d.get(i).scale(alpha.get(i, a)));

						rarray[a] = rarray[a].minus(
								p.get(i).scale(alpha.get(i, a)));

					}

					if (Utils.mag(rarray[a]) < 1E-6) {//todo change this if you want
						rarray[a] = null;
					}
					else {
//						System.out.println("convergence test: " + mag
//								(rarray[a]));
					}
				}
			}

			solver = new SimpleMatrix(solver.numRows(),
					solver.numRows());

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					double[] arrrhs = new double[solver.numRows()];

					for (int i = 0; i < arrrhs.length; i++) {
						arrrhs[i] = -rarray[a].transpose()
								.mult(p.get(i)).get(0, 0);

					}
					rhsvec.setColumn(a, 0, arrrhs);
				}
			}

			for (int i = 0; i < solver.numRows(); i++) {
				for (int j = 0; j < solver.numRows(); j++) {
					solver.set(i, j,
							d.get(j).transpose().mult(p.get(i))
									.get(0, 0));
				}
			}

			SimpleMatrix beta = solver.solve(rhsvec);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					dirs[a] = rarray[a];

					for (int i = 0; i < beta.numRows(); i++) {
						dirs[a] = dirs[a].plus(
								d.get(i).scale(beta.get(i, a)));
					}
				}
			}
		}

		SimpleMatrix[] densityMatrixDerivs =
				new SimpleMatrix[fockderivstatic.length];

		for (int a = 0; a < fockderivstatic.length; a++) {
			SimpleMatrix densityMatrixDeriv = new SimpleMatrix
					(soln.orbitals.length, soln.orbitals.length);

			for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
				for (int v = u; v < densityMatrixDeriv.numCols(); v++) {
					double sum = 0;
					int count = 0;
					for (int i = 0; i < NOcc; i++) {
						for (int j = 0; j < NVirt; j++) {
							sum -= 2 * (soln.C.get(i, u) *
									soln.C.get(j + NOcc, v) +
									soln.C.get(j + NOcc, u) *
											soln.C.get(i, v)) *
									xarray[a].get(count, 0);
							count++;
						}
					}

					densityMatrixDeriv.set(u, v, sum);
					densityMatrixDeriv.set(v, u, sum);
				}
			}

			densityMatrixDerivs[a] = densityMatrixDeriv;
		}

		return densityMatrixDerivs;
	}

	public static SimpleMatrix[] getxarrayPople(SolutionR soln,
												SimpleMatrix[] fockderivstatic) {
		int NOcc = (int) (soln.nElectrons / 2.0);
		int NVirt = soln.orbitals.length - NOcc;
		int length = fockderivstatic.length;
		int nonv = NOcc * NVirt;

		if (nonv == 0) {
			SimpleMatrix[] xarray = new SimpleMatrix[length];

			for (int i = 0; i < xarray.length; i++) {
				xarray[i] = new SimpleMatrix(0, 0);
			}

			return xarray;
		}

		// array initialization
		SimpleMatrix[] xarray = new SimpleMatrix[length];
		SimpleMatrix[] barray = new SimpleMatrix[length];
		SimpleMatrix[] parray = new SimpleMatrix[length];
		SimpleMatrix[] Farray = new SimpleMatrix[length];
		SimpleMatrix[] rarray = new SimpleMatrix[length];

		// configure preconditioners
		double[] Darr = new double[nonv];
		double[] Dinvarr = new double[nonv];

		int counter = 0;
		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = -soln.E.get(i) + soln.E.get(NOcc + j);

				Darr[counter] = Math.pow(e, -0.5);
				Dinvarr[counter] = Math.pow(e, 0.5);

				counter++;
			}
		}

		// convert AO to MO basis
		SimpleMatrix F = new SimpleMatrix(nonv, length);
		for (int a = 0; a < length; a++) {
			SimpleMatrix f = new SimpleMatrix(nonv, 1);

			int count = 0;

			for (int i = 0; i < NOcc; i++) { // kappa
				for (int j = 0; j < NVirt; j++) { // i
					double element = 0;

					for (int u = 0; u < soln.orbitals.length; u++) {
						for (int v = 0; v < soln.orbitals.length; v++) {
							element += soln.C.get(i, u) *
									soln.C.get(j + NOcc, v) *
									fockderivstatic[a].get(u, v);
						}
					}

					element /= soln.E.get(j + NOcc) - soln.E.get(i);

					f.set(count, 0, element);

					count++;
				}
			}

			CommonOps_DDRM.multRows(Darr, f.getDDRM());
			barray[a] = f;
			Farray[a] = barray[a].copy();
			F.setColumn(a, 0, barray[a].getDDRM().data);
		}

		// main loop
		int[] iterable = new int[length];

		// 0: B, 1: Bt, 2: Bn, 3: P, 4: BmP
		List<SimpleMatrix[]> prevs = new ArrayList<>();
		List<Double> dots = new ArrayList<>();

		while (Utils.numIterable(iterable) > 0) {
			// orthogonalize barray
			for (int i = 0; i < barray.length; i++) {
				for (int j = 0; j < i; j++) {
					barray[i].plusi(barray[i].dot(barray[j]) /
							barray[j].dot(barray[j]), barray[j].negativei());
				}
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix[] prev = new SimpleMatrix[5];
				prev[0] = barray[i]; // original barray object here
				prev[1] = barray[i].transpose();
				prev[2] = barray[i].negative();
				dots.add(barray[i].dot(barray[i]));

				// parray[i] stays the same object throughout
				SimpleMatrix bc = barray[i].copy();
				CommonOps_DDRM.multRows(Dinvarr, bc.getDDRM());
				SimpleMatrix crv = computeResponseVectorsPople(bc, soln);
				CommonOps_DDRM.multRows(Darr, crv.getDDRM());
				parray[i] = crv;

				prev[3] = parray[i];
				prev[4] = barray[i].minus(parray[i]);

				prevs.add(prev);
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix newb = parray[i].copy();

				// orthogonalize against all previous Bs
				for (int j = 0; j < prevs.size(); j++) {
					SimpleMatrix[] prev = prevs.get(j);
					SimpleMatrix transpose = prev[1];
					double num = transpose.mult(parray[i]).get(0) /
							dots.get(j);

					newb.plusi(num, prev[2]);
				}

				barray[i] = newb; // new barray object created
			}

			SimpleMatrix Bt = new SimpleMatrix(prevs.size(), nonv);
			SimpleMatrix P = new SimpleMatrix(nonv, prevs.size());

			for (int i = 0; i < prevs.size(); i++) {
				Bt.setRow(i, 0, prevs.get(i)[0].getDDRM().data);
				P.setColumn(i, 0, prevs.get(i)[4].getDDRM().data);
			}

			SimpleMatrix rhs = Bt.mult(F);
			SimpleMatrix lhs = Bt.mult(P);
			// alpha dimensions are prevBs x length
			SimpleMatrix alpha = lhs.solve(rhs);

			// reset r and x array
			for (int a = 0; a < length; a++) {
				rarray[a] = new SimpleMatrix(nonv, 1);
				xarray[a] = new SimpleMatrix(nonv, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {
					// B with tilde
					rarray[j].plusi(alpha.get(i, j), prevs.get(i)[4]);
					xarray[j].plusi(alpha.get(i, j), prevs.get(i)[0]);
				}
			}

			for (int j = 0; j < alpha.numCols(); j++) {
				// B0 is Farray, no tilde
				rarray[j].minusi(Farray[j]);
				CommonOps_DDRM.multRows(Dinvarr, xarray[j].getDDRM());

				double rMag = Utils.mag(rarray[j]);
				if (rMag < 1E-7) {
					iterable[j] = 1;
				}
				else if (Double.isNaN(rMag)) {
					soln.getRm().getLogger()
							.warn("Pople algorithm fails; reverting to " +
									"Thiel algorithm (don't panic)...");
					throw new SingularMatrixException();
				}
				else {
					iterable[j] = 0;
				}
			}
		}

		return xarray;
	}

	public static SimpleMatrix[] densityDerivPople(SolutionR soln,
												   SimpleMatrix[] fockderivstatic) {
		int NOcc = (int) (soln.nElectrons / 2.0);
		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix[] xarray = getxarrayPople(soln, fockderivstatic);

		SimpleMatrix[] densityMatrixDerivs =
				new SimpleMatrix[fockderivstatic.length];

		for (int a = 0; a < fockderivstatic.length; a++) {
			SimpleMatrix densityMatrixDeriv =
					new SimpleMatrix(soln.orbitals.length,
							soln.orbitals.length);

			for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
				for (int v = u; v < densityMatrixDeriv.numCols(); v++) {
					double sum = 0;
					int count = 0;
					for (int i = 0; i < NOcc; i++) {
						for (int j = 0; j < NVirt; j++) {
							sum -= 2 * (soln.C.get(i, u) *
									soln.C.get(j + NOcc, v) +
									soln.C.get(j + NOcc, u) *
											soln.C.get(i, v)) *
									xarray[a].get(count, 0);
							count++;
						}
					}

					densityMatrixDeriv.set(u, v, sum);
					densityMatrixDeriv.set(v, u, sum);
				}
			}

			densityMatrixDerivs[a] = densityMatrixDeriv;
		}

		return densityMatrixDerivs;
	}


	private static SimpleMatrix computeResponseVectorsThiel(SimpleMatrix x,
															SolutionR soln) {

		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix densityMatrixDeriv =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
			for (int v = u; v < densityMatrixDeriv.numCols(); v++) {
				double sum = 0;
				int count = 0;
				for (int i = 0; i < NOcc; i++) {
					for (int j = 0; j < NVirt; j++) {
						sum -= 2 * (soln.C.get(i, u) * soln.C.get(j + NOcc,
								v) +
								soln.C.get(j + NOcc, u) * soln.C.get(i, v)) *
								x.get(count, 0);
						count++;
					}
				}

				densityMatrixDeriv.set(u, v, sum);
				densityMatrixDeriv.set(v, u, sum);
			}
		}


		SimpleMatrix responsematrix =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = j; k < soln.orbitals.length; k++) {
				double val = 0;
				if (j == k) {

					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							val += densityMatrixDeriv.get(l, l) *
									integralArray[integralcount];
							integralcount++;
						}
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] ==
											soln.atomOfOrb[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += densityMatrixDeriv.get(j, k) *
							integralArray[integralcount];
					integralcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] ==
											soln.atomOfOrb[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.orbsOfAtom[soln.atomOfOrb[k]]) {
								if (m > -1) {
									val += densityMatrixDeriv.get(l, m) *
											integralArray[integralcount];
									integralcount++;
								}
							}
						}
					}
				}

				responsematrix.set(j, k, val);
				responsematrix.set(k, j, val);
			}
		}

		SimpleMatrix R = new SimpleMatrix(NOcc * NVirt, 1);

		int count1 = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {

				double element = 0;

				for (int u = 0; u < soln.orbitals.length; u++) {
					for (int v = 0; v < soln.orbitals.length; v++) {
						element += soln.C.get(i, u) * soln.C.get(j + NOcc, v) *
								responsematrix.get(u, v);
					}
				}


				R.set(count1, 0, element);

				count1++;
			}
		}


		SimpleMatrix p = new SimpleMatrix(NOcc * NVirt, 1);

		int counter = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				p.set(counter, 0, -R.get(counter, 0) +
						(soln.E.get(j + NOcc) - soln.E.get(i)) *
								x.get(counter));
				counter++;
			}
		}


		return p;
	}

	public static SimpleMatrix computeResponseVectorsPople(SimpleMatrix x,
														   SolutionR soln) {
		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix densityMatrixDeriv =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
			for (int v = u; v < densityMatrixDeriv.numCols(); v++) {
				double sum = 0;
				int count = 0;
				for (int i = 0; i < NOcc; i++) {
					for (int j = 0; j < NVirt; j++) {
						sum -= 2 * (soln.C.get(i, u) * soln.C.get(j + NOcc,
								v) +
								soln.C.get(j + NOcc, u) * soln.C.get(i, v)) *
								x.get(count, 0);
						count++;
					}
				}

				densityMatrixDeriv.set(u, v, sum);
				densityMatrixDeriv.set(v, u, sum);
			}
		}


		SimpleMatrix responsematrix =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = j; k < soln.orbitals.length; k++) {
				double val = 0;
				if (j == k) {

					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							val += densityMatrixDeriv.get(l, l) *
									integralArray[integralcount];
							integralcount++;
						}
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] ==
											soln.atomOfOrb[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += densityMatrixDeriv.get(j, k) *
							integralArray[integralcount];
					integralcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] ==
											soln.atomOfOrb[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.orbsOfAtom[soln.atomOfOrb[k]]) {
								if (m > -1) {
									val += densityMatrixDeriv.get(l, m) *
											integralArray[integralcount];
									integralcount++;
								}
							}
						}
					}
				}

				responsematrix.set(j, k, val);
				responsematrix.set(k, j, val);
			}
		}

		SimpleMatrix R = new SimpleMatrix(NOcc * NVirt, 1);

		int count1 = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {

				double element = 0;

				for (int u = 0; u < soln.orbitals.length; u++) {
					for (int v = 0; v < soln.orbitals.length; v++) {
						element += soln.C.get(i, u) * soln.C.get(j + NOcc, v) *
								responsematrix.get(u, v);
					}
				}

				element = element / (soln.E.get(j + NOcc) - soln.E.get(i));


				R.set(count1, 0, element);

				count1++;
			}
		}


		return R;
	}
}
