package nddoparam.param;

import nddoparam.*;
import org.jblas.DoubleMatrix;
import org.jblas.Solve;
import org.jblas.exceptions.LapackException;
import scf.GTO;
import scf.LCGTO;
import scf.Utils;

import java.lang.reflect.Constructor;
import java.util.ArrayList;

public class ParamDerivative {

	private static double p1Deriv(NDDOAtom a, int type) {

		double D1deriv = D1Deriv(a, type);

		if (a.getAtomProperties().getZ() == 1) {
			return 0;
		}

		return -a.p1 * a.p1 * a.D1 /
				(a.p1 * a.p1 * a.p1 -
						Math.pow(a.D1 * a.D1 + a.p1 * a.p1, 1.5)) * D1deriv;
	}


	private static double p2Deriv(NDDOAtom a, int type) {

		if (type == 0) {
			return 0;
		}

		double D2deriv = D2Deriv(a, type);

		double F1 = 2 * a.D2 * (Math.pow(a.D2 * a.D2 + a.p2 * a.p2, -1.5) -
				Math.pow(2 * a.D2 * a.D2 + a.p2 * a.p2, -1.5));

		double F2 = a.p2 * (2 * Math.pow(a.D2 * a.D2 + a.p2 * a.p2, -1.5) -
				Math.pow(2 * a.D2 * a.D2 + a.p2 * a.p2, -1.5)) -
				1 / (a.p2 * a.p2);

		return -F1 / F2 * D2deriv;

	}

	private static double D1Deriv(NDDOAtom a, int type) {

		double zetas = a.getParams().getZetas();

		double zetap = a.getParams().getZetap();

		double zeta = 0;

		switch (type) {

			case 0:
				zeta = zetap;
				break;
			case 1:
				zeta = zetas;
				break;
			default:
				zeta = 0;
		}

		return (2 * a.getAtomProperties().getPeriod() + 1) / Math.sqrt(3) *
				(4 * zeta * (0.5 + a.getAtomProperties().getPeriod()) *
						Math.pow(4 * zetas * zetap,
								-0.5 + a.getAtomProperties().getPeriod()) /
						Math.pow(zetas + zetap,
								2 + 2 * a.getAtomProperties().getPeriod())
						- (2 + 2 * a.getAtomProperties().getPeriod()) *
						Math.pow(4 * zetas * zetap,
								0.5 + a.getAtomProperties().getPeriod()) /
						Math.pow(zetas + zetap,
								3 + 2 * a.getAtomProperties().getPeriod()));


	}

	private static double D2Deriv(NDDOAtom a, int type) {
		if (type == 0) {
			return 0;
		}

		return -1 / a.getParams().getZetap() * a.D2;
	}

	private static double qq(double p01, double p11, double p21, double D11,
							 double D21,
							 double p02, double p12, double p22, double D12,
							 double D22,
							 double R, int num, double D1deriv, double D2deriv,
							 double p1deriv, double p2deriv) {
		return 0;
	}

	private static double quz(double p01, double p11, double p21, double D11,
							  double D21,
							  double p02, double p12, double p22, double D12,
							  double D22,
							  double R, int num, double D1deriv,
							  double D2deriv,
							  double p1deriv, double p2deriv) {
		double a01 = p01 + p12;

		if (num == 0) {
			return 0;
		}
		if (num == 1) {
			return -0.5 * ((R + D12) * D1deriv + a01 * p1deriv) *
					Math.pow((R + D12) * (R + D12) + a01 * a01, -1.5)
					+ 0.5 * ((D12 - R) * D1deriv + a01 * p1deriv) *
					Math.pow((R - D12) * (R - D12) + a01 * a01, -1.5);
		}
		if (num == 2) {
			return -0.5 * ((R + D11) * D1deriv + a01 * p1deriv) *
					Math.pow((R + D11) * (R + D11) + a01 * a01, -1.5)
					+ 0.5 * ((D11 - R) * D1deriv + a01 * p1deriv) *
					Math.pow((R - D11) * (R - D11) + a01 * a01, -1.5);
		}

		return 0;
	}

	private static double qQpipi(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R, int num,
								 double D1deriv, double D2deriv,
								 double p1deriv,
								 double p2deriv) {
		double a02 = p01 + p22;

		if (num == 0) {
			return 0;
		}
		if (num == 1) {
			return -0.5 * (4 * D22 * D2deriv + a02 * p2deriv) *
					Math.pow(R * R + 4 * D22 * D22 + a02 * a02, -1.5)
					+ 0.5 * a02 * p2deriv * Math.pow(R * R + a02 * a02, -1.5);
		}
		if (num == 2) {
			return -0.5 * (4 * D21 * D2deriv + a02 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a02 * a02, -1.5)
					+ 0.5 * a02 * p2deriv * Math.pow(R * R + a02 * a02, -1.5);
		}

		return 0;

	}

	private static double qQzz(double p01, double p11, double p21, double D11,
							   double D21,
							   double p02, double p12, double p22, double D12,
							   double D22,
							   double R, int num, double D1deriv,
							   double D2deriv,
							   double p1deriv, double p2deriv) {
		double a02 = p01 + p22;

		if (num == 0) {
			return 0;
		}
		if (num == 1) {
			return -0.25 * (2 * (R + 2 * D22) * D2deriv + a02 * p2deriv) *
					Math.pow((R + 2 * D22) * (R + 2 * D22) + a02 * a02, -1.5)
					- 0.25 * (2 * (2 * D22 - R) * D2deriv + a02 * p2deriv) *
					Math.pow((R - 2 * D22) * (R - 2 * D22) + a02 * a02, -1.5)
					+ 0.5 * a02 * p2deriv * Math.pow(R * R + a02 * a02, -1.5);
		}

		if (num == 2) {
			return -0.25 * (2 * (R + 2 * D21) * D2deriv + a02 * p2deriv) *
					Math.pow((R + 2 * D21) * (R + 2 * D21) + a02 * a02, -1.5)
					- 0.25 * (2 * (2 * D21 - R) * D2deriv + a02 * p2deriv) *
					Math.pow((R - 2 * D21) * (R - 2 * D21) + a02 * a02, -1.5)
					+ 0.5 * a02 * p2deriv * Math.pow(R * R + a02 * a02, -1.5);
		}

		return 0;

	}

	private static double upiupi(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R, int num,
								 double D1deriv, double D2deriv,
								 double p1deriv,
								 double p2deriv) {
		double a11 = p11 + p12;
		if (num == 0) {
			return -0.5 * ((D11 - D12) * D1deriv + a11 * p1deriv) *
					Math.pow(R * R + (D11 - D12) * (D11 - D12) + a11 * a11,
							-1.5)
					+ 0.5 * ((D11 + D12) * D1deriv + a11 * p1deriv) *
					Math.pow(R * R + (D11 + D12) * (D11 + D12) + a11 * a11,
							-1.5);
		}
		if (num == 1) {
			return -0.5 * ((D12 - D11) * D1deriv + a11 * p1deriv) *
					Math.pow(R * R + (D11 - D12) * (D11 - D12) + a11 * a11,
							-1.5)
					+ 0.5 * ((D11 + D12) * D1deriv + a11 * p1deriv) *
					Math.pow(R * R + (D11 + D12) * (D11 + D12) + a11 * a11,
							-1.5);
		}
		if (num == 2) {
			return -0.5 * (2 * a11 * p1deriv) *
					Math.pow(R * R + a11 * a11, -1.5)
					+ 0.5 * (4 * D11 * D1deriv + 2 * a11 * p1deriv) *
					Math.pow(R * R + 4 * D11 * D11 + a11 * a11, -1.5);
		}

		return 0;

	}

	private static double uzuz(double p01, double p11, double p21, double D11,
							   double D21,
							   double p02, double p12, double p22, double D12,
							   double D22,
							   double R, int num, double D1deriv,
							   double D2deriv,
							   double p1deriv, double p2deriv) {
		double a11 = p11 + p12;

		if (num == 0) {
			return -0.25 * ((R + D11 - D12) * D1deriv + a11 * p1deriv) *
					Math.pow((R + D11 - D12) * (R + D11 - D12) + a11 * a11,
							-1.5)
					+ 0.25 * ((R + D11 + D12) * D1deriv + a11 * p1deriv) *
					Math.pow((R + D11 + D12) * (R + D11 + D12) + a11 * a11,
							-1.5)
					+ 0.25 * ((D11 + D12 - R) * D1deriv + a11 * p1deriv) *
					Math.pow((R - D11 - D12) * (R - D11 - D12) + a11 * a11,
							-1.5)
					- 0.25 * ((D11 - D12 - R) * D1deriv + a11 * p1deriv) *
					Math.pow((R - D11 + D12) * (R - D11 + D12) + a11 * a11,
							-1.5);
		}
		if (num == 1) {
			return -0.25 * ((D12 - D11 - R) * D1deriv + a11 * p1deriv) *
					Math.pow((R + D11 - D12) * (R + D11 - D12) + a11 * a11,
							-1.5)
					+ 0.25 * ((R + D11 + D12) * D1deriv + a11 * p1deriv) *
					Math.pow((R + D11 + D12) * (R + D11 + D12) + a11 * a11,
							-1.5)
					+ 0.25 * ((D11 + D12 - R) * D1deriv + a11 * p1deriv) *
					Math.pow((R - D11 - D12) * (R - D11 - D12) + a11 * a11,
							-1.5)
					- 0.25 * ((D12 + R - D11) * D1deriv + a11 * p1deriv) *
					Math.pow((R - D11 + D12) * (R - D11 + D12) + a11 * a11,
							-1.5);
		}

		if (num == 2) {
			return -a11 * p1deriv * Math.pow(R * R + a11 * a11, -1.5)
					+ 0.25 * (2 * (R + 2 * D11) * D1deriv + 2 * a11 * p1deriv) *
					Math.pow((R + 2 * D11) * (R + 2 * D11) + a11 * a11, -1.5)
					+ 0.25 * (2 * (2 * D11 - R) * D1deriv + 2 * a11 * p1deriv) *
					Math.pow((R - 2 * D11) * (R - 2 * D11) + a11 * a11, -1.5);
		}

		return 0;

	}

	private static double upiQpiz(double p01, double p11, double p21,
								  double D11,
								  double D21, double p02, double p12,
								  double p22,
								  double D12, double D22, double R, int num,
								  double D1deriv, double D2deriv,
								  double p1deriv,
								  double p2deriv) {
		double a12 = p11 + p22;

		if (num == 0) {
			return 0.25 * ((D11 - D22) * D1deriv + a12 * p1deriv) * Math.pow(
					(R - D22) * (R - D22) + (D11 - D22) * (D11 - D22) +
							a12 * a12, -1.5)
					- 0.25 * ((D11 + D22) * D1deriv + a12 * p1deriv) * Math.pow(
					(R - D22) * (R - D22) + (D11 + D22) * (D11 + D22) +
							a12 * a12, -1.5)
					- 0.25 * ((D11 - D22) * D1deriv + a12 * p1deriv) * Math.pow(
					(R + D22) * (R + D22) + (D11 - D22) * (D11 - D22) +
							a12 * a12, -1.5)
					+ 0.25 * ((D11 + D22) * D1deriv + a12 * p1deriv) * Math.pow(
					(R + D22) * (R + D22) + (D11 + D22) * (D11 + D22) +
							a12 * a12, -1.5);
		}
		if (num == 1) {
			return 0.25 * ((2 * D22 - R - D11) * D2deriv + a12 * p2deriv) *
					Math.pow(
							(R - D22) * (R - D22) + (D11 - D22) * (D11 - D22) +
									a12 * a12, -1.5)
					- 0.25 * ((2 * D22 + D11 - R) * D2deriv + a12 * p2deriv) *
					Math.pow(
							(R - D22) * (R - D22) + (D11 + D22) * (D11 + D22) +
									a12 * a12, -1.5)
					- 0.25 * ((2 * D22 + R - D11) * D2deriv + a12 * p2deriv) *
					Math.pow(
							(R + D22) * (R + D22) + (D11 - D22) * (D11 - D22) +
									a12 * a12, -1.5)
					+ 0.25 * ((2 * D22 + D11 + R) * D2deriv + a12 * p2deriv) *
					Math.pow(
							(R + D22) * (R + D22) + (D11 + D22) * (D11 + D22) +
									a12 * a12, -1.5);
		}

		if (num == 2) {
			return 0.25 *
					((D21 - R) * D2deriv + (D11 - D21) * (D1deriv - D2deriv) +
							a12 * (p1deriv + p2deriv)) * Math.pow(
					(R - D21) * (R - D21) + (D11 - D21) * (D11 - D21) +
							a12 * a12, -1.5)
					- 0.25 *
					((D21 - R) * D2deriv + (D11 + D21) * (D1deriv + D2deriv) +
							a12 * (p1deriv + p2deriv)) * Math.pow(
					(R - D21) * (R - D21) + (D11 + D21) * (D11 + D21) +
							a12 * a12, -1.5)
					- 0.25 *
					((D21 + R) * D2deriv + (D11 - D21) * (D1deriv - D2deriv) +
							a12 * (p1deriv + p2deriv)) * Math.pow(
					(R + D21) * (R + D21) + (D11 - D21) * (D11 - D21) +
							a12 * a12, -1.5)
					+ 0.25 *
					((D21 + R) * D2deriv + (D11 + D21) * (D1deriv + D2deriv) +
							a12 * (p1deriv + p2deriv)) * Math.pow(
					(R + D21) * (R + D21) + (D11 + D21) * (D11 + D21) +
							a12 * a12, -1.5);
		}

		return 0;

	}

	private static double uzQpipi(double p01, double p11, double p21,
								  double D11,
								  double D21, double p02, double p12,
								  double p22,
								  double D12, double D22, double R, int num,
								  double D1deriv, double D2deriv,
								  double p1deriv,
								  double p2deriv) {
		double a12 = p11 + p22;

		if (num == 0) {
			return 0.25 * ((R + D11) * D1deriv + a12 * p1deriv) *
					Math.pow((R + D11) * (R + D11) + 4 * D22 * D22 + a12 * a12,
							-1.5)
					- 0.25 * ((D11 - R) * D1deriv + a12 * p1deriv) *
					Math.pow((R - D11) * (R - D11) + 4 * D22 * D22 + a12 * a12,
							-1.5)
					- 0.25 * ((R + D11) * D1deriv + a12 * p1deriv) *
					Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
					+ 0.25 * ((D11 - R) * D1deriv + a12 * p1deriv) *
					Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		}
		if (num == 1) {
			return 0.25 * (4 * D22 * D2deriv + a12 * p2deriv) *
					Math.pow((R + D11) * (R + D11) + 4 * D22 * D22 + a12 * a12,
							-1.5)
					- 0.25 * (4 * D22 * D2deriv + a12 * p2deriv) *
					Math.pow((R - D11) * (R - D11) + 4 * D22 * D22 + a12 * a12,
							-1.5)
					- 0.25 * (a12 * p2deriv) *
					Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
					+ 0.25 * (a12 * p2deriv) *
					Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		}

		if (num == 2) {
			return 0.25 * ((R + D11) * D1deriv + 4 * D21 * D2deriv +
					a12 * (p1deriv + p2deriv)) *
					Math.pow((R + D11) * (R + D11) + 4 * D21 * D21 + a12 * a12,
							-1.5)
					- 0.25 * ((D11 - R) * D1deriv + 4 * D21 * D2deriv +
					a12 * (p1deriv + p2deriv)) *
					Math.pow((R - D11) * (R - D11) + 4 * D21 * D21 + a12 * a12,
							-1.5)
					- 0.25 * ((R + D11) * D1deriv + a12 * (p1deriv + p2deriv)) *
					Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
					+ 0.25 * ((D11 - R) * D1deriv + a12 * (p1deriv + p2deriv)) *
					Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		}

		return 0;

	}

	private static double uzQzz(double p01, double p11, double p21, double D11,
								double D21, double p02, double p12, double p22,
								double D12, double D22, double R, int num,
								double D1deriv,
								double D2deriv, double p1deriv,
								double p2deriv) {
		double a12 = p11 + p22;

		if (num == 0) {
			return +0.125 * ((R + D11 - 2 * D22) * D1deriv + a12 * p1deriv) *
					Math.pow((R + D11 - 2 * D22) * (R + D11 - 2 * D22) +
							a12 * a12, -1.5)
					- 0.125 * ((D11 + 2 * D22 - R) * D1deriv + a12 * p1deriv) *
					Math.pow((R - D11 - 2 * D22) * (R - D11 - 2 * D22) +
							a12 * a12, -1.5)
					+ 0.125 * ((R + D11 + 2 * D22) * D1deriv + a12 * p1deriv) *
					Math.pow((R + D11 + 2 * D22) * (R + D11 + 2 * D22) +
							a12 * a12, -1.5)
					- 0.125 * ((D11 - 2 * D22 - R) * D1deriv + a12 * p1deriv) *
					Math.pow((R - D11 + 2 * D22) * (R - D11 + 2 * D22) +
							a12 * a12, -1.5)
					- 0.25 * ((R + D11) * D1deriv + a12 * p1deriv) *
					Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
					+ 0.25 * ((D11 - R) * D1deriv + a12 * p1deriv) *
					Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		}
		if (num == 1) {
			return +0.125 *
					(2 * (2 * D22 - D11 - R) * D2deriv + a12 * p2deriv) *
					Math.pow((R + D11 - 2 * D22) * (R + D11 - 2 * D22) +
							a12 * a12, -1.5)
					- 0.125 *
					(2 * (D11 + 2 * D22 - R) * D2deriv + a12 * p2deriv) *
					Math.pow((R - D11 - 2 * D22) * (R - D11 - 2 * D22) +
							a12 * a12, -1.5)
					+ 0.125 *
					(2 * (R + D11 + 2 * D22) * D2deriv + a12 * p2deriv) *
					Math.pow((R + D11 + 2 * D22) * (R + D11 + 2 * D22) +
							a12 * a12, -1.5)
					- 0.125 *
					(2 * (2 * D22 + R - D11) * D2deriv + a12 * p2deriv) *
					Math.pow((R - D11 + 2 * D22) * (R - D11 + 2 * D22) +
							a12 * a12, -1.5)
					- 0.25 * (a12 * p2deriv) *
					Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
					+ 0.25 * (a12 * p2deriv) *
					Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		}

		if (num == 2) {
			return +0.125 * ((R + D11 - 2 * D21) * (D1deriv - 2 * D2deriv) +
					a12 * (p1deriv + p2deriv)) *
					Math.pow((R + D11 - 2 * D21) * (R + D11 - 2 * D21) +
							a12 * a12, -1.5)
					- 0.125 * ((D11 + 2 * D21 - R) * (D1deriv + 2 * D2deriv) +
					a12 * (p1deriv + p2deriv)) *
					Math.pow((R - D11 - 2 * D21) * (R - D11 - 2 * D21) +
							a12 * a12, -1.5)
					+ 0.125 * ((R + D11 + 2 * D21) * (D1deriv + 2 * D2deriv) +
					a12 * (p1deriv + p2deriv)) *
					Math.pow((R + D11 + 2 * D21) * (R + D11 + 2 * D21) +
							a12 * a12, -1.5)
					- 0.125 * ((R - D11 + 2 * D21) * (2 * D2deriv - D1deriv) +
					a12 * (p1deriv + p2deriv)) *
					Math.pow((R - D11 + 2 * D21) * (R - D11 + 2 * D21) +
							a12 * a12, -1.5)
					- 0.25 * ((R + D11) * D1deriv + a12 * (p1deriv + p2deriv)) *
					Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
					+ 0.25 * ((D11 - R) * D1deriv + a12 * (p1deriv + p2deriv)) *
					Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		}

		return 0;

	}

	private static double QpipiQpipi(double p01, double p11, double p21,
									 double D11,
									 double D21, double p02, double p12,
									 double p22,
									 double D12, double D22, double R, int num,
									 double D1deriv, double D2deriv,
									 double p1deriv,
									 double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) {
			return -0.125 * (4 * (D21 - D22) * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * (D21 - D22) * (D21 - D22) + a22 * a22,
							-1.5)
					- 0.125 * (4 * (D21 + D22) * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * (D21 + D22) * (D21 + D22) + a22 * a22,
							-1.5)
					+ 0.25 * (4 * D21 * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
					+ 0.25 * (a22 * p2deriv) *
					Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5)
					-
					0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		}
		if (num == 1) {
			return -0.125 * (4 * (D22 - D21) * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * (D21 - D22) * (D21 - D22) + a22 * a22,
							-1.5)
					- 0.125 * (4 * (D21 + D22) * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * (D21 + D22) * (D21 + D22) + a22 * a22,
							-1.5)
					+ 0.25 * (a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
					+ 0.25 * (4 * D22 * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5)
					-
					0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		}

		if (num == 2) {
			return -0.75 * a22 * p2deriv * Math.pow(R * R + a22 * a22, -1.5)
					- 0.125 * (16 * D21 * D2deriv + 2 * a22 * p2deriv) *
					Math.pow(R * R + 16 * D21 * D21 + a22 * a22, -1.5)
					+ 0.5 * (4 * D21 * D2deriv + 2 * a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5);
		}

		return 0;
	}

	private static double QxxQyy(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R, int num,
								 double D1deriv, double D2deriv,
								 double p1deriv,
								 double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) {
			return -0.25 * (4 * D21 * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22,
							-1.5)
					+ 0.25 * (4 * D21 * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
					+ 0.25 * (a22 * p2deriv) *
					Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5)
					-
					0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		}
		if (num == 1) {
			return -0.25 * (4 * D22 * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22,
							-1.5)
					+ 0.25 * (a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
					+ 0.25 * (4 * D22 * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5)
					-
					0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		}

		if (num == 2) {
			return -0.25 * (8 * D22 * D2deriv + 2 * a22 * p2deriv) *
					Math.pow(R * R + 8 * D21 * D21 + a22 * a22, -1.5)
					+ 0.5 * (4 * D22 * D2deriv + 2 * a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
					- 0.25 * (2 * a22 * p2deriv) *
					Math.pow(R * R + a22 * a22, -1.5);
		}

		return 0;

	}

	private static double QpipiQzz(double p01, double p11, double p21,
								   double D11,
								   double D21, double p02, double p12,
								   double p22,
								   double D12, double D22, double R, int num,
								   double D1deriv, double D2deriv,
								   double p1deriv,
								   double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) {
			return -0.125 * (4 * D21 * D2deriv + a22 * p2deriv) *
					Math.pow((R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 +
									a22 * a22,
							-1.5)
					- 0.125 * (4 * D21 * D2deriv + a22 * p2deriv) *
					Math.pow((R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 +
									a22 * a22,
							-1.5)
					+ 0.125 * (a22 * p2deriv) *
					Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5)
					+ 0.125 * (a22 * p2deriv) *
					Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5)
					+ 0.25 * (4 * D21 * D2deriv + a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
					-
					0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		}
		if (num == 1) {
			return -0.125 * (2 * (2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 +
									a22 * a22,
							-1.5)
					- 0.125 * (2 * (2 * D22 + R) * D2deriv + a22 * p2deriv) *
					Math.pow((R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 +
									a22 * a22,
							-1.5)
					+ 0.125 * (2 * (2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5)
					+ 0.125 * (2 * (2 * D22 + R) * D2deriv + a22 * p2deriv) *
					Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5)
					+ 0.25 * (a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
					-
					0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		}

		if (num == 2) {
			return -0.125 *
					((2 * (2 * D21 - R) + 4 * D21) * D2deriv +
							2 * a22 * p2deriv) *
					Math.pow((R - 2 * D21) * (R - 2 * D21) + 4 * D21 * D21 +
									a22 * a22,
							-1.5)
					- 0.125 *
					((2 * (2 * D21 + R) + 4 * D21) * D2deriv +
							2 * a22 * p2deriv) *
					Math.pow((R + 2 * D21) * (R + 2 * D21) + 4 * D21 * D21 +
									a22 * a22,
							-1.5)
					+
					0.125 * (2 * (2 * D21 - R) * D2deriv + 2 * a22 * p2deriv) *
							Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22,
									-1.5)
					+
					0.125 * (2 * (2 * D21 + R) * D2deriv + 2 * a22 * p2deriv) *
							Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22,
									-1.5)
					+ 0.25 * (4 * D21 * D2deriv + 2 * a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
					- 0.25 * (2 * a22 * p2deriv) *
					Math.pow(R * R + a22 * a22, -1.5);
		}

		return 0;

	}

	private static double QzzQzz(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R, int num,
								 double D1deriv, double D2deriv,
								 double p1deriv,
								 double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) {

			return -0.0625 *
					(2 * (R + 2 * D21 - 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) +
									a22 * a22,
							-1.5)
					- 0.0625 *
					(2 * (R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) +
									a22 * a22,
							-1.5)
					- 0.0625 *
					(2 * (2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) +
									a22 * a22,
							-1.5)
					- 0.0625 *
					(2 * (2 * D21 - R - 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) +
									a22 * a22,
							-1.5)
					+ 0.125 * (2 * (R + 2 * D21) * D2deriv + a22 * p2deriv) *
					Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5)
					+ 0.125 * (2 * (2 * D21 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5)
					+ 0.125 * (a22 * p2deriv) *
					Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5)
					+ 0.125 * (a22 * p2deriv) *
					Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5)
					-
					0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		}
		if (num == 1) {


			return -0.0625 *
					(2 * (2 * D22 - R - 2 * D21) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) +
									a22 * a22,
							-1.5)
					- 0.0625 *
					(2 * (R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) +
									a22 * a22,
							-1.5)
					- 0.0625 *
					(2 * (2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) +
									a22 * a22,
							-1.5)
					- 0.0625 *
					(2 * (R - 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) +
									a22 * a22,
							-1.5)
					+ 0.125 * (a22 * p2deriv) *
					Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5)
					+ 0.125 * (a22 * p2deriv) *
					Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5)
					+ 0.125 * (2 * (R + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5)
					+ 0.125 * (2 * (2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5)
					-
					0.25 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		}

		if (num == 2) {
			return -0.0625 * (4 * (R + 4 * D21) * D2deriv + 2 * a22 * p2deriv) *
					Math.pow((R + 4 * D21) * (R + 4 * D21) + a22 * a22, -1.5)
					-
					0.0625 * (4 * (4 * D21 - R) * D2deriv + 2 * a22 * p2deriv) *
							Math.pow((R - 4 * D21) * (R - 4 * D21) + a22 * a22,
									-1.5)
					+ 0.25 * (2 * (R + 2 * D21) * D2deriv + 2 * a22 * p2deriv) *
					Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5)
					+ 0.25 * (2 * (2 * D21 - R) * D2deriv + 2 * a22 * p2deriv) *
					Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5)
					-
					0.75 * (a22 * p2deriv) * Math.pow(R * R + a22 * a22, -1.5);
		}

		return 0;
	}

	private static double QpizQpiz(double p01, double p11, double p21,
								   double D11,
								   double D21, double p02, double p12,
								   double p22,
								   double D12, double D22, double R, int num,
								   double D1deriv, double D2deriv,
								   double p1deriv,
								   double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) {
			return -0.125 *
					((R + 2 * D21 - 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 - D22) * (R + D21 - D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					+ 0.125 * ((R + 2 * D21) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R + D21 - D22) * (R + D21 - D22) +
									(D21 + D22) * (D21 + D22) +
									a22 * a22, -1.5)
					+ 0.125 * ((R + 2 * D21) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R + D21 + D22) * (R + D21 + D22) +
									(D21 - D22) * (D21 - D22) +
									a22 * a22, -1.5)
					- 0.125 *
					((R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 + D22) * (R + D21 + D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					+ 0.125 * ((2 * D21 - R) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R - D21 - D22) * (R - D21 - D22) +
									(D21 - D22) * (D21 - D22) +
									a22 * a22, -1.5)
					- 0.125 *
					((2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 - D22) * (R - D21 - D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					- 0.125 *
					((2 * D21 - 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 + D22) * (R - D21 + D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					+ 0.125 * ((2 * D21 - R) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R - D21 + D22) * (R - D21 + D22) +
									(D21 + D22) * (D21 + D22) +
									a22 * a22, -1.5);
		}
		if (num == 1) {
			return -0.125 *
					((2 * D22 - 2 * D21 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 - D22) * (R + D21 - D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					+ 0.125 * ((2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R + D21 - D22) * (R + D21 - D22) +
									(D21 + D22) * (D21 + D22) +
									a22 * a22, -1.5)
					+ 0.125 * ((R + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R + D21 + D22) * (R + D21 + D22) +
									(D21 - D22) * (D21 - D22) +
									a22 * a22, -1.5)
					- 0.125 *
					((R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 + D22) * (R + D21 + D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					+ 0.125 * ((2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R - D21 - D22) * (R - D21 - D22) +
									(D21 - D22) * (D21 - D22) +
									a22 * a22, -1.5)
					- 0.125 *
					((2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 - D22) * (R - D21 - D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					- 0.125 *
					((2 * D22 + R - 2 * D21) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 + D22) * (R - D21 + D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					+ 0.125 * ((R + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow(
							(R - D21 + D22) * (R - D21 + D22) +
									(D21 + D22) * (D21 + D22) +
									a22 * a22, -1.5);
		}

		if (num == 2) {
			return -0.5 * a22 * p2deriv * Math.pow(R * R + a22 * a22, -1.5)
					+ 0.25 * (4 * D21 * D2deriv + 2 * a22 * p2deriv) *
					Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
					+
					0.125 * (2 * (R + 2 * D21) * D2deriv + 2 * a22 * p2deriv) *
							Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22,
									-1.5)
					+
					0.125 * (2 * (2 * D21 - R) * D2deriv + 2 * a22 * p2deriv) *
							Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22,
									-1.5)
					- 0.125 *
					((2 * (2 * D21 + R) + 4 * D21) * D2deriv +
							2 * a22 * p2deriv) *
					Math.pow((R + 2 * D21) * (R + 2 * D21) + 4 * D21 * D21 +
									a22 * a22,
							-1.5)
					- 0.125 *
					((2 * (2 * D21 - R) + 4 * D21) * D2deriv +
							2 * a22 * p2deriv) *
					Math.pow((R - 2 * D21) * (R - 2 * D21) + 4 * D21 * D21 +
									a22 * a22,
							-1.5);


		}

		return 0;
	}

	private static double ssssderiv(double p01, double p11, double p21,
									double D11,
									double D21, double p02, double p12,
									double p22,
									double D12, double D22, double R, int num,
									double D1deriv, double D2deriv,
									double p1deriv,
									double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv);
	}

	private static double ssppippideriv(double p01, double p11, double p21,
										double D11,
										double D21, double p02, double p12,
										double p22,
										double D12, double D22, double R,
										int num,
										double D1deriv, double D2deriv,
										double p1deriv,
										double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num,
						D1deriv,
						D2deriv, p1deriv, p2deriv);
	}

	private static double sspzpzderiv(double p01, double p11, double p21,
									  double D11,
									  double D21, double p02, double p12,
									  double p22,
									  double D12, double D22, double R,
									  int num,
									  double D1deriv, double D2deriv,
									  double p1deriv,
									  double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
						D1deriv,
						D2deriv, p1deriv, p2deriv);
	}

	private static double ppippissderiv(double p01, double p11, double p21,
										double D11,
										double D21, double p02, double p12,
										double p22,
										double D12, double D22, double R,
										int num,
										double D1deriv, double D2deriv,
										double p1deriv,
										double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pzpzssderiv(double p01, double p11, double p21,
									  double D11,
									  double D21, double p02, double p12,
									  double p22,
									  double D12, double D22, double R,
									  int num,
									  double D1deriv, double D2deriv,
									  double p1deriv,
									  double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv,
						D2deriv, p1deriv, p2deriv);
	}

	private static double ppippippippideriv(double p01, double p11, double p21,
											double D11, double D21, double p02,
											double p12, double p22, double D12,
											double D22, double R, int num,
											double D1deriv,
											double D2deriv, double p1deriv,
											double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num,
						D1deriv,
						D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv, D2deriv, p1deriv, p2deriv) +
				QpipiQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num,
						D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pxpxpypyderiv(double p01, double p11, double p21,
										double D11,
										double D21, double p02, double p12,
										double p22,
										double D12, double D22, double R,
										int num,
										double D1deriv, double D2deriv,
										double p1deriv,
										double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num,
						D1deriv,
						D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv, D2deriv, p1deriv, p2deriv) +
				QxxQyy(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num,
						D1deriv,
						D2deriv, p1deriv, p2deriv);
	}

	private static double ppippipzpzderiv(double p01, double p11, double p21,
										  double D11,
										  double D21, double p02, double p12,
										  double p22,
										  double D12, double D22, double R,
										  int num,
										  double D1deriv, double D2deriv,
										  double p1deriv,
										  double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
						D1deriv,
						D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv, D2deriv, p1deriv, p2deriv) +
				QpipiQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num,
						D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pzpzppippideriv(double p01, double p11, double p21,
										  double D11,
										  double D21, double p02, double p12,
										  double p22,
										  double D12, double D22, double R,
										  int num,
										  double D1deriv, double D2deriv,
										  double p1deriv,
										  double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num,
						D1deriv,
						D2deriv, p1deriv, p2deriv) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv,
						D2deriv, p1deriv, p2deriv) +
				QpipiQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pzpzpzpzderiv(double p01, double p11, double p21,
										double D11,
										double D21, double p02, double p12,
										double p22,
										double D12, double D22, double R,
										int num,
										double D1deriv, double D2deriv,
										double p1deriv,
										double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
						D1deriv,
						D2deriv, p1deriv, p2deriv) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv,
						D2deriv, p1deriv, p2deriv) +
				QzzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num,
						D1deriv,
						D2deriv, p1deriv, p2deriv);
	}

	private static double spzssderiv(double p01, double p11, double p21,
									 double D11,
									 double D21, double p02, double p12,
									 double p22,
									 double D12, double D22, double R, int num,
									 double D1deriv, double D2deriv,
									 double p1deriv,
									 double p2deriv) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
				f(num),
				D1deriv,
				D2deriv, p1deriv, p2deriv);
	}

	private static double spzppippideriv(double p01, double p11, double p21,
										 double D11,
										 double D21, double p02, double p12,
										 double p22,
										 double D12, double D22, double R,
										 int num,
										 double D1deriv, double D2deriv,
										 double p1deriv,
										 double p2deriv) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
				f(num),
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				uzQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num,
						D1deriv,
						D2deriv, p1deriv, p2deriv);
	}

	private static double spzpzpzderiv(double p01, double p11, double p21,
									   double D11,
									   double D21, double p02, double p12,
									   double p22,
									   double D12, double D22, double R,
									   int num,
									   double D1deriv, double D2deriv,
									   double p1deriv,
									   double p2deriv) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
				f(num),
				D1deriv,
				D2deriv, p1deriv, p2deriv) +
				uzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
						D1deriv,
						D2deriv, p1deriv, p2deriv);
	}

	private static double ssspzderiv(double p01, double p11, double p21,
									 double D11,
									 double D21, double p02, double p12,
									 double p22,
									 double D12, double D22, double R, int num,
									 double D1deriv, double D2deriv,
									 double p1deriv,
									 double p2deriv) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv);
	}

	private static double ppippispzderiv(double p01, double p11, double p21,
										 double D11,
										 double D21, double p02, double p12,
										 double p22,
										 double D12, double D22, double R,
										 int num,
										 double D1deriv, double D2deriv,
										 double p1deriv,
										 double p2deriv) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) -
				uzQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pzpzspzderiv(double p01, double p11, double p21,
									   double D11,
									   double D21, double p02, double p12,
									   double p22,
									   double D12, double D22, double R,
									   int num,
									   double D1deriv, double D2deriv,
									   double p1deriv,
									   double p2deriv) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv) -
				uzQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num),
						D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double sppisppideriv(double p01, double p11, double p21,
										double D11,
										double D21, double p02, double p12,
										double p22,
										double D12, double D22, double R,
										int num,
										double D1deriv, double D2deriv,
										double p1deriv,
										double p2deriv) {
		return upiupi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv);
	}

	private static double spzspzderiv(double p01, double p11, double p21,
									  double D11,
									  double D21, double p02, double p12,
									  double p22,
									  double D12, double D22, double R,
									  int num,
									  double D1deriv, double D2deriv,
									  double p1deriv,
									  double p2deriv) {
		return uzuz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv,
				D2deriv, p1deriv, p2deriv);
	}

	private static double sppippipzderiv(double p01, double p11, double p21,
										 double D11,
										 double D21, double p02, double p12,
										 double p22,
										 double D12, double D22, double R,
										 int num,
										 double D1deriv, double D2deriv,
										 double p1deriv,
										 double p2deriv) {
		return upiQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
				num,
				D1deriv,
				D2deriv, p1deriv, p2deriv);
	}

	private static double ppipzsppideriv(double p01, double p11, double p21,
										 double D11,
										 double D21, double p02, double p12,
										 double p22,
										 double D12, double D22, double R,
										 int num,
										 double D1deriv, double D2deriv,
										 double p1deriv,
										 double p2deriv) {
		return -upiQpiz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
				f(num),
				D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double ppipzppipzderiv(double p01, double p11, double p21,
										  double D11,
										  double D21, double p02, double p12,
										  double p22,
										  double D12, double D22, double R,
										  int num,
										  double D1deriv, double D2deriv,
										  double p1deriv,
										  double p2deriv) {
		return QpizQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
				num,
				D1deriv,
				D2deriv, p1deriv, p2deriv);
	}

	private static double pxpypxpyderiv(double p01, double p11, double p21,
										double D11,
										double D21, double p02, double p12,
										double p22,
										double D12, double D22, double R,
										int num,
										double D1deriv, double D2deriv,
										double p1deriv,
										double p2deriv) {
		return 0.5 *
				(ppippippippideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12,
						D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv) -
						pxpxpypyderiv(p01, p11, p21, D11, D21, p02, p12, p22,
								D12, D22
								, R,
								num, D1deriv, D2deriv, p1deriv, p2deriv));
	}

	private static int f(int num) {

		if (num == 0 || num == 1) {
			return 1 - num;
		}
		else if (num == 2) {
			return num;
		}

		return -1;
	}


	protected static double LocalTwoCenterERIderiv(NDDO6G a, NDDO6G b,
												   NDDO6G c,
												   NDDO6G d,
												   double D1deriv,
												   double D2deriv,
												   double p1deriv,
												   double p2deriv,
												   int num, int type) {


		double[] A = a.getCoords();
		double[] C = c.getCoords();

		double R = GTO.R(A, C);
		//(??|??)
		switch (a.getL()) {

			case 0://(s?|??)

				switch (b.getL()) {

					case 0: //(ss|??)

						switch (c.getL()) {

							case 0: //(ss|s?);

								switch (d.getL()) {

									case 0://(ss|ss)
										return ssssderiv(a.p0, a.p1, a.p2,
												a.D1,
												a.D2,
												c.p0, c.p1, c.p2, c.D1, c.D2
												, R,
												num,
												D1deriv, D2deriv, p1deriv,
												p2deriv);

									case 1:
										if (d.getk() == 1) {//(ss|spz)
											return ssspzderiv(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);
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
											return ssspzderiv(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);

										case 1:
											if (d.getk() == 1) {//(ss|pzpz)
												return sspzpzderiv(a.p0, a.p1,
														a.p2,
														a.D1,
														a.D2, c.p0, c.p1, c.p2,
														c.D1,
														c.D2, R, num, D1deriv,
														D2deriv,
														p1deriv, p2deriv);
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
										return ssppippideriv(a.p0, a.p1, a.p2,
												a.D1,
												a.D2,
												c.p0, c.p1, c.p2, c.D1, c.D2
												, R,
												num,
												D1deriv, D2deriv, p1deriv,
												p2deriv);
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
											return spzssderiv(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);

										case 1:
											if (d.getk() == 1) {//(spz|spz)
												return spzspzderiv(a.p0, a.p1,
														a.p2,
														a.D1,
														a.D2, c.p0, c.p1, c.p2,
														c.D1,
														c.D2, R, num, D1deriv,
														D2deriv,
														p1deriv, p2deriv);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(spz|pz?)

										switch (d.getL()) {

											case 0://(spz|pzs)
												return spzspzderiv(a.p0, a.p1,
														a.p2,
														a.D1,
														a.D2, c.p0, c.p1, c.p2,
														c.D1,
														c.D2, R, num, D1deriv,
														D2deriv,
														p1deriv, p2deriv);

											case 1:
												if (d.getk() == 1) {//(spz
													// |pzpz)
													return spzpzpzderiv(a.p0,
															a.p1, a.p2,
															a.D1, a.D2, c.p0,
															c.p1, c.p2,
															c.D1, c.D2, R, num,
															D1deriv,
															D2deriv, p1deriv,
															p2deriv);
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
											return spzppippideriv(a.p0, a.p1,
													a.p2, a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);
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
										return sppisppideriv(a.p0, a.p1, a.p2,
												a.D1,
												a.D2,
												c.p0, c.p1, c.p2, c.D1, c.D2
												, R,
												num,
												D1deriv, D2deriv, p1deriv,
												p2deriv);
									}
									else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == b.geti() &&
												d.getj() == b.getj() &&
												d.getk() == 0) {//(sppi|pzppi)
											return sppippipzderiv(a.p0, a.p1,
													a.p2, a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);
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
													return sppisppideriv(a.p0,
															a.p1,
															a.p2,
															a.D1, a.D2, c.p0,
															c.p1, c.p2,
															c.D1, c.D2, R, num,
															D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												case 1:
													if (d.getk() == 1) {
														return sppippipzderiv(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R,
																num, D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
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
											return spzssderiv(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);

										case 1:
											if (d.getk() == 1) {//(pzs|spz)
												return spzspzderiv(a.p0, a.p1,
														a.p2,
														a.D1,
														a.D2, c.p0, c.p1, c.p2,
														c.D1,
														c.D2, R, num, D1deriv,
														D2deriv,
														p1deriv, p2deriv);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(pzs|pz?)

										switch (d.getL()) {

											case 0://(pzs|pzs)
												return spzspzderiv(a.p0, a.p1,
														a.p2,
														a.D1,
														a.D2, c.p0, c.p1, c.p2,
														c.D1,
														c.D2, R, num, D1deriv,
														D2deriv,
														p1deriv, p2deriv);

											case 1:
												if (d.getk() == 1) {//(pzs
													// |pzpz)
													return spzpzpzderiv(a.p0,
															a.p1, a.p2,
															a.D1, a.D2, c.p0,
															c.p1, c.p2,
															c.D1, c.D2, R, num,
															D1deriv,
															D2deriv, p1deriv,
															p2deriv);
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
											return spzppippideriv(a.p0, a.p1,
													a.p2, a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);
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
												return pzpzssderiv(a.p0, a.p1,
														a.p2,
														a.D1,
														a.D2, c.p0, c.p1, c.p2,
														c.D1,
														c.D2, R, num, D1deriv,
														D2deriv,
														p1deriv, p2deriv);

											case 1:
												if (d.getk() == 1) {//(pzpz
													// |spz)
													return pzpzspzderiv(a.p0,
															a.p1, a.p2,
															a.D1, a.D2, c.p0,
															c.p1, c.p2,
															c.D1, c.D2, R, num,
															D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												}
												else {
													return 0;
												}
										}

									case 1:
										if (c.getk() == 1) {//(pzpz|pz?)

											switch (d.getL()) {

												case 0://(pzpz|pzs)
													return pzpzspzderiv(a.p0,
															a.p1, a.p2,
															a.D1, a.D2, c.p0,
															c.p1, c.p2,
															c.D1, c.D2, R, num,
															D1deriv,
															D2deriv, p1deriv,
															p2deriv);

												case 1:
													if (d.getk() ==
															1) {//(pzpz|pzpz)
														return pzpzpzpzderiv(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R,
																num, D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
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
												return pzpzppippideriv(a.p0,
														a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2,
														c.D1, c.D2, R, num,
														D1deriv,
														D2deriv, p1deriv,
														p2deriv);
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
											return ppipzsppideriv(a.p0, a.p1,
													a.p2, a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);
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
												return ppipzppipzderiv(a.p0,
														a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2,
														c.D1, c.D2, R, num,
														D1deriv,
														D2deriv, p1deriv,
														p2deriv);
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
														return ppipzsppideriv(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R,
																num, D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzderiv(
																	a.p0,
																	a.p1, a.p2,
																	a.D1,
																	a.D2, c.p0,
																	c.p1,
																	c.p2, c.D1,
																	c.D2, R,
																	num,
																	D1deriv,
																	D2deriv,
																	p1deriv,
																	p2deriv);
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
										return sppisppideriv(a.p0, a.p1, a.p2,
												a.D1,
												a.D2,
												c.p0, c.p1, c.p2, c.D1, c.D2
												, R,
												num,
												D1deriv, D2deriv, p1deriv,
												p2deriv);
									}
									else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == a.geti() &&
												d.getj() == a.getj() &&
												d.getk() == 0) {//(ppis|pzppi)
											return sppippipzderiv(a.p0, a.p1,
													a.p2, a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);
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
													return sppisppideriv(a.p0,
															a.p1,
															a.p2,
															a.D1, a.D2, c.p0,
															c.p1, c.p2,
															c.D1, c.D2, R, num,
															D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												case 1:
													if (d.getk() == 1) {
														return sppippipzderiv(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R,
																num, D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
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
											return ppipzsppideriv(a.p0, a.p1,
													a.p2, a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2
													, R,
													num, D1deriv, D2deriv,
													p1deriv,
													p2deriv);
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
												return ppipzppipzderiv(a.p0,
														a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2,
														c.D1, c.D2, R, num,
														D1deriv,
														D2deriv, p1deriv,
														p2deriv);
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
														return ppipzsppideriv(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R,
																num, D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzderiv(
																	a.p0,
																	a.p1, a.p2,
																	a.D1,
																	a.D2, c.p0,
																	c.p1,
																	c.p2, c.D1,
																	c.D2, R,
																	num,
																	D1deriv,
																	D2deriv,
																	p1deriv,
																	p2deriv);
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
													return ppippissderiv(a.p0,
															a.p1,
															a.p2,
															a.D1, a.D2, c.p0,
															c.p1, c.p2,
															c.D1, c.D2, R, num,
															D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												}
												else {
													return 0;
												}
											case 1:
												if (d.getk() == 1 &&
														a.geti() == b.geti() &&
														a.getj() == b.getj()) {
													return ppippispzderiv(a.p0,
															a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R,
															num,
															D1deriv, D2deriv,
															p1deriv,
															p2deriv);
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
														return ppippispzderiv(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R,
																num, D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													}
													else {
														return 0;
													}

												case 1:
													if (d.getk() == 1 &&
															a.geti() ==
																	b.geti() &&
															a.getj() ==
																	b.getj()) {//(ppippi
														// |pzpz)
														return ppippipzpzderiv(
																a.p0,
																a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R,
																num, D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
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
														return ppippippippideriv(
																a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													}
													else {
														return pxpxpypyderiv(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R,
																num, D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
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
													return pxpypxpyderiv(a.p0,
															a.p1,
															a.p2,
															a.D1, a.D2, c.p0,
															c.p1, c.p2,
															c.D1, c.D2, R, num,
															D1deriv,
															D2deriv, p1deriv,
															p2deriv);
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


	public static double getGderiv(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d,
								   int num,
								   int type) {
		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double sum2 = 0;

		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		double p1deriv = 0;

		double p2deriv = 0;

		double D1deriv = 0;

		double D2deriv = 0;

		if (num == 0 || num == 2) {

			p1deriv = p1Deriv(a.getAtom(), type);
			p2deriv = p2Deriv(a.getAtom(), type);
			D1deriv = D1Deriv(a.getAtom(), type);
			D2deriv = D2Deriv(a.getAtom(), type);
		}
		else if (num == 1) {

			p1deriv = p1Deriv(c.getAtom(), type);
			p2deriv = p2Deriv(c.getAtom(), type);
			D1deriv = D1Deriv(c.getAtom(), type);
			D2deriv = D2Deriv(c.getAtom(), type);
		}

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {

						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] !=
								0) {
							sum2 += coeffA[i] * coeffB[j] * coeffC[k] *
									coeffD[l] *
									LocalTwoCenterERIderiv(A[i], B[j], C[k],
											D[l],
											D1deriv, D2deriv, p1deriv, p2deriv,
											num,
											type) * 27.21;
						}


					}
				}
			}
		}


		return sum2;
	}


	public static double getGderivfinite(NDDO6G a, NDDO6G b, NDDO6G c,
										 NDDO6G d,
										 int num,
										 int type) {

		int aindex = index(a);

		int bindex = index(b);

		int cindex = index(c);

		int dindex = index(d);

		double initial = NDDO6G.getG(a, b, c, d);

		NDDOAtom A = a.getAtom();

		NDDOAtom C = c.getAtom();

		if (num == 0 || num == 2) {

			try {
				NDDOParams params = A.getParams().clone();
				params.modifyParam(5 + type, Utils.LAMBDA);

				Class<? extends NDDOAtom> cl = A.getClass();
				Constructor ctor =
						cl.getDeclaredConstructor(cl,
								A.getParams().getClass());
				ctor.setAccessible(true);

				A = (NDDOAtom) ctor.newInstance(A, params);
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(0);
			}
		}
		if (num == 1 || num == 2) {
			try {
				NDDOParams params = C.getParams().clone();
				params.modifyParam(5 + type, Utils.LAMBDA);

				Class<? extends NDDOAtom> cl = C.getClass();
				Constructor ctor =
						cl.getDeclaredConstructor(cl,
								C.getParams().getClass());
				ctor.setAccessible(true);

				C = (NDDOAtom) ctor.newInstance(C, params);
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(0);
			}
		}

		double finalval =
				NDDO6G.getG(A.getOrbitals()[aindex], A.getOrbitals()[bindex],
						C.getOrbitals()[cindex], C.getOrbitals()[dindex]);

		return (finalval - initial) / Utils.LAMBDA;


	}

	public static double getSderivfinite(NDDO6G a, NDDO6G b, int num,
										 int type) {

		int aindex = index(a);

		int bindex = index(b);

		double initial = NDDO6G.getS(a, b);

		NDDOAtom A = a.getAtom();

		NDDOAtom B = b.getAtom();

		if (num == 0 || num == 2) {

			try {
				NDDOParams params = A.getParams().clone();
				params.modifyParam(5 + type, Utils.LAMBDA);

				Class<? extends NDDOAtom> cl = A.getClass();
				Constructor ctor =
						cl.getDeclaredConstructor(cl,
								A.getParams().getClass());
				ctor.setAccessible(true);

				A = (NDDOAtom) ctor.newInstance(A, params);
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(0);
			}
		}
		if (num == 1 || num == 2) {
			try {
				NDDOParams params = B.getParams().clone();
				params.modifyParam(5 + type, Utils.LAMBDA);

				Class<? extends NDDOAtom> cl = B.getClass();
				Constructor ctor =
						cl.getDeclaredConstructor(cl,
								B.getParams().getClass());
				ctor.setAccessible(true);

				B = (NDDOAtom) ctor.newInstance(B, params);
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(0);
			}
		}

		double finalval =
				NDDO6G.getS(A.getOrbitals()[aindex], B.getOrbitals()[bindex]);

		return (finalval - initial) / Utils.LAMBDA;


	}

	public static double crfderivfinite(NDDOAtom A, NDDOAtom B, int num) {

		double initial = A.crf(B);

		if (num == 0 || num == 2) {

			try {
				NDDOParams params = A.getParams().clone();
				params.modifyParam(0, Utils.LAMBDA);

				Class<? extends NDDOAtom> cl = A.getClass();
				Constructor ctor =
						cl.getDeclaredConstructor(cl,
								A.getParams().getClass());
				ctor.setAccessible(true);

				A = (NDDOAtom) ctor.newInstance(A, params);
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(0);
			}
		}
		if (num == 1 || num == 2) {
			try {
				NDDOParams params = B.getParams().clone();
				params.modifyParam(0, Utils.LAMBDA);

				Class<? extends NDDOAtom> cl = B.getClass();
				Constructor ctor =
						cl.getDeclaredConstructor(cl,
								B.getParams().getClass());
				ctor.setAccessible(true);

				B = (NDDOAtom) ctor.newInstance(B, params);
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(0);
			}
		}

		double finalval = A.crf(B);

		return (finalval - initial) / Utils.LAMBDA;


	}

	public static double HFDeriv(SolutionR soln, int Z, int paramnum) {

		if (paramnum == 0) {
			return alphaHfderiv(soln, Z);
		}
		if (paramnum <= 2) {
			return betaHfderiv(soln, Z, paramnum - 1);
		}
		if (paramnum <= 4) {
			return uxxHfderiv(soln, Z, paramnum - 3);
		}

		if (paramnum <= 6) {
			return zetaHfderiv(soln, Z, paramnum - 5);
		}
		if (paramnum == 7) {
			return eisolHfderiv(soln.atoms, Z);
		}

		System.err.println("oh no! This isn't MNDO!");
		System.exit(0);
		return 0;

	}

	public static double[] MNDOHfderivs(SolutionR soln, int Z) {

		double[] derivs = new double[8];

		derivs[0] = alphaHfderiv(soln, Z);
		derivs[1] = betaHfderiv(soln, Z, 0);
		derivs[3] = uxxHfderiv(soln, Z, 0);
		derivs[5] = zetaHfderiv(soln, Z, 0);
		derivs[7] = eisolHfderiv(soln.atoms, Z);

		if (Z != 1) {
			derivs[2] = betaHfderiv(soln, Z, 1);
			derivs[4] = uxxHfderiv(soln, Z, 1);
			derivs[6] = zetaHfderiv(soln, Z, 1);

			return derivs;
		}
		else {
			return new double[]{derivs[0], derivs[1], derivs[3], derivs[5],
					derivs[7]};
		}
	}


	public static DoubleMatrix[][] MNDOStaticMatrixDeriv(SolutionR soln,
														 int Z,
														 int firstParamIndex) {
		NDDOAtom[] atoms = soln.atoms;
		DoubleMatrix[] HDerivs = new DoubleMatrix[8];
		DoubleMatrix[] FDerivs = new DoubleMatrix[8];

		if (firstParamIndex <= 1) HDerivs[1] = betafockderivstatic(soln, Z, 0);
		if (firstParamIndex <= 1) FDerivs[1] = HDerivs[1].dup();
		if (firstParamIndex <= 3) HDerivs[3] = uxxfockderivstatic(soln, Z, 0);
		if (firstParamIndex <= 3) FDerivs[3] = HDerivs[3].dup();
		if (firstParamIndex <= 5)
			HDerivs[5] = zetaHderivstatic(atoms, soln, Z, 0);
		if (firstParamIndex <= 5)
			FDerivs[5] =
					HDerivs[5].dup().add(zetaGderivstatic(atoms, soln, Z, 0));

		if (Z != 1) {
			if (firstParamIndex <= 2)
				HDerivs[2] = betafockderivstatic(soln, Z, 1);
			if (firstParamIndex <= 2) FDerivs[2] = HDerivs[2].dup();
			if (firstParamIndex <= 4)
				HDerivs[4] = uxxfockderivstatic(soln, Z, 1);
			if (firstParamIndex <= 4) FDerivs[4] = HDerivs[4].dup();
			if (firstParamIndex <= 6)
				HDerivs[6] = zetaHderivstatic(atoms, soln, Z, 1);
			if (firstParamIndex <= 6)
				FDerivs[6] = HDerivs[6].dup()
						.add(zetaGderivstatic(atoms, soln, Z, 1));
		}
		return new DoubleMatrix[][]{HDerivs, FDerivs};
	}

	public static double MNDOHFDeriv(SolutionR soln, DoubleMatrix Hderiv,
									 DoubleMatrix Fderiv) {

		double e = 0;

		DoubleMatrix densitymatrix = soln.densityMatrix();

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = 0; k < soln.orbitals.length; k++) {
				e += 0.5 * densitymatrix.get(j, k) *
						(Hderiv.get(j, k) + Fderiv.get(j, k));
			}
		}

		return e / 4.3363E-2;
	}

	private static double zetaHfderiv(SolutionR soln, int Z, int type) {

		DoubleMatrix densitymatrix = soln.densityMatrix();


		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.index;

		int[][] missingIndex = soln.missingIndex;

		int[] atomNumber = soln.atomNumber;

		int[] atomicnumbers = soln.atomicNumbers;

		DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] == atomNumber[k]) {
					double Huv = 0;

					for (int an = 0; an < soln.atoms.length; an++) {
						if (atomNumber[j] != an) {
							Huv += soln.atoms[an]
									.VParamDeriv(orbitals[j], orbitals[k],
											getNum(atomicnumbers[an],
													atomicnumbers[atomNumber[j]],
													Z), type);
						}
					}
					H.put(j, k, Huv);
					H.put(k, j, Huv);
				}
				else { // case 3
					double Huk = NDDO6G.betaparamderiv(orbitals[j],
							orbitals[k],
							getNum(atomicnumbers[atomNumber[j]],
									atomicnumbers[atomNumber[k]], Z), type);
					H.put(j, k, Huk);
					H.put(k, j, Huk);
				}
			}
		}

		DoubleMatrix G = DoubleMatrix.zeros(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;

				if (atomNumber[j] == atomNumber[k]) {
					for (int l : missingIndex[atomNumber[j]]) {
						if (l > -1) {
							for (int m : missingIndex[atomNumber[j]]) {
								if (m > -1) {
									if (atomNumber[l] == atomNumber[m]) {
										sum += soln.densityMatrix().get(l, m) *
												(ParamDerivative
														.getGderiv(orbitals[j],
																orbitals[k],
																orbitals[l],
																orbitals[m],
																getNum(atomicnumbers[atomNumber[j]],
																		atomicnumbers[atomNumber[l]],
																		Z),
																type));

									}
								}
							}
						}
					}
				}
				else {
					for (int l : index[atomNumber[j]]) {
						if (l > -1) {
							for (int m : index[atomNumber[k]]) {
								if (m > -1) {
									sum += soln.densityMatrix().get(l, m) *
											(-0.5 *
													ParamDerivative
															.getGderiv(
																	orbitals[j],
																	orbitals[l],
																	orbitals[k],
																	orbitals[m],
																	getNum(atomicnumbers[atomNumber[j]],
																			atomicnumbers[atomNumber[k]],
																			Z),
																	type));
								}
							}
						}
					}
				}

				G.put(j, k, sum);
				G.put(k, j, sum);
			}
		}

		DoubleMatrix F = H.dup().add(G);

		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * densitymatrix.get(j, k) *
						(H.get(j, k) + F.get(j, k));
			}
		}

		return e / 4.3363E-2;

	}

	public static DoubleMatrix zetafockderivstatic(NDDOAtom[] atoms,
												   SolutionR soln, int Z,
												   int type) {


		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.index;

		int[][] missingIndex = soln.missingIndex;

		int[] atomNumber = soln.atomNumber;

		int[] atomicnumbers = soln.atomicNumbers;

		DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] == atomNumber[k]) {
					double Huv = 0;

					for (int an = 0; an < atoms.length; an++) {
						if (atomNumber[j] != an) {
							Huv += atoms[an]
									.VParamDeriv(orbitals[j], orbitals[k],
											getNum(atomicnumbers[an],
													atomicnumbers[atomNumber[j]],
													Z), type);
						}
					}
					H.put(j, k, Huv);
					H.put(k, j, Huv);
				}
				else { // case 3
					double Huk = NDDO6G.betaparamderiv(orbitals[j],
							orbitals[k],
							getNum(atomicnumbers[atomNumber[j]],
									atomicnumbers[atomNumber[k]], Z), type);
					H.put(j, k, Huk);
					H.put(k, j, Huk);
				}
			}
		}

		DoubleMatrix G = DoubleMatrix.zeros(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;

				if (atomNumber[j] == atomNumber[k]) {
					for (int l : missingIndex[atomNumber[j]]) {
						if (l > -1) {
							for (int m : missingIndex[atomNumber[j]]) {
								if (m > -1) {
									if (atomNumber[l] == atomNumber[m]) {
										sum += soln.densityMatrix().get(l, m) *
												(ParamDerivative
														.getGderiv(orbitals[j],
																orbitals[k],
																orbitals[l],
																orbitals[m],
																getNum(atomicnumbers[atomNumber[j]],
																		atomicnumbers[atomNumber[l]],
																		Z),
																type));

									}
								}
							}
						}
					}
				}
				else {
					for (int l : index[atomNumber[j]]) {
						if (l > -1) {
							for (int m : index[atomNumber[k]]) {
								if (m > -1) {
									sum += soln.densityMatrix().get(l, m) *
											(-0.5 *
													ParamDerivative
															.getGderiv(
																	orbitals[j],
																	orbitals[l],
																	orbitals[k],
																	orbitals[m],
																	getNum(atomicnumbers[atomNumber[j]],
																			atomicnumbers[atomNumber[k]],
																			Z),
																	type));
								}
							}
						}
					}
				}

				G.put(j, k, sum);
				G.put(k, j, sum);
			}
		}

		DoubleMatrix F = H.dup().add(G);

		return F;
	}

	public static DoubleMatrix zetaHderivstatic(NDDOAtom[] atoms,
												SolutionR soln, int Z,
												int type) {


		NDDO6G[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomNumber;

		int[] atomicnumbers = soln.atomicNumbers;

		DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] == atomNumber[k]) {
					double Huv = 0;

					for (int an = 0; an < atoms.length; an++) {
						if (atomNumber[j] != an) {
							Huv += atoms[an]
									.VParamDeriv(orbitals[j], orbitals[k],
											getNum(atomicnumbers[an],
													atomicnumbers[atomNumber[j]],
													Z), type);
						}
					}
					H.put(j, k, Huv);
					H.put(k, j, Huv);
				}
				else { // case 3
					double Huk = NDDO6G.betaparamderiv(orbitals[j],
							orbitals[k],
							getNum(atomicnumbers[atomNumber[j]],
									atomicnumbers[atomNumber[k]], Z), type);
					H.put(j, k, Huk);
					H.put(k, j, Huk);
				}
			}
		}

		return H;
	}

	public static DoubleMatrix zetaGderivstatic(NDDOAtom[] atoms,
												SolutionR soln, int Z,
												int type) {


		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.index;

		int[][] missingIndex = soln.missingIndex;

		int[] atomNumber = soln.atomNumber;

		int[] atomicnumbers = soln.atomicNumbers;

		DoubleMatrix G = DoubleMatrix.zeros(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;

				if (atomNumber[j] == atomNumber[k]) {
					for (int l : missingIndex[atomNumber[j]]) {
						if (l > -1) {
							for (int m : missingIndex[atomNumber[j]]) {
								if (m > -1) {
									if (atomNumber[l] == atomNumber[m]) {
										sum += soln.densityMatrix().get(l, m) *
												(ParamDerivative
														.getGderiv(orbitals[j],
																orbitals[k],
																orbitals[l],
																orbitals[m],
																getNum(atomicnumbers[atomNumber[j]],
																		atomicnumbers[atomNumber[l]],
																		Z),
																type));

									}
								}
							}
						}
					}
				}
				else {
					for (int l : index[atomNumber[j]]) {
						if (l > -1) {
							for (int m : index[atomNumber[k]]) {
								if (m > -1) {
									sum += soln.densityMatrix().get(l, m) *
											(-0.5 *
													ParamDerivative
															.getGderiv(
																	orbitals[j],
																	orbitals[l],
																	orbitals[k],
																	orbitals[m],
																	getNum(atomicnumbers[atomNumber[j]],
																			atomicnumbers[atomNumber[k]],
																			Z),
																	type));
								}
							}
						}
					}
				}

				G.put(j, k, sum);
				G.put(k, j, sum);
			}
		}


		return G;
	}


	private static double uxxHfderiv(SolutionR soln, int Z, int type) {

		DoubleMatrix densitymatrix = soln.densityMatrix();


		NDDO6G[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomNumber;

		int[] atomicnumbers = soln.atomicNumbers;

		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {

			if (atomicnumbers[atomNumber[j]] == Z &&
					orbitals[j].getL() == type) {
				e += densitymatrix.get(j, j);
			}
		}

		return e / 4.3363E-2;

	}

	public static DoubleMatrix uxxfockderivstatic(SolutionR soln, int Z,
												  int type) {

		DoubleMatrix F =
				DoubleMatrix.zeros(soln.orbitals.length, soln.orbitals.length);

		for (int j = 0; j < soln.orbitals.length; j++) {
			if (soln.atomicNumbers[soln.atomNumber[j]] == Z &&
					soln.orbitals[j].getL() == type) {
				F.put(j, j, 1);
			}
		}

		return F;

	}

	private static double betaHfderiv(SolutionR soln, int Z, int type) {

		DoubleMatrix densitymatrix = soln.densityMatrix();


		NDDO6G[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomNumber;

		int[] atomicnumbers = soln.atomicNumbers;


		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] != atomNumber[k]) {

					double H = betaderiv(atomicnumbers[atomNumber[j]],
							atomicnumbers[atomNumber[k]], Z,
							orbitals[j].getL(),
							orbitals[k].getL(), type);

					if (H != 0) {
						H *= LCGTO.getS(orbitals[j], orbitals[k]);
					}
					e += 2 * densitymatrix.get(j, k) * H;
				}
			}
		}

		return e / 4.3363E-2;

	}

	public static DoubleMatrix betafockderivstatic(SolutionR soln, int Z,
												   int type) {

		DoubleMatrix F =
				DoubleMatrix.zeros(soln.orbitals.length, soln.orbitals.length);

		NDDO6G[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomNumber;

		int[] atomicnumbers = soln.atomicNumbers;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] != atomNumber[k]) {

					double H = betaderiv(atomicnumbers[atomNumber[j]],
							atomicnumbers[atomNumber[k]], Z,
							orbitals[j].getL(),
							orbitals[k].getL(), type);

					if (H != 0) {
						H *= LCGTO.getS(orbitals[j], orbitals[k]);
					}
					F.put(j, k, H);
					F.put(k, j, H);
				}
			}
		}

		return F;

	}


	private static double eisolHfderiv(NDDOAtom[] atoms, int Z) {

		int counter = 0;

		for (NDDOAtom a : atoms) {
			if (a.getAtomProperties().getZ() == Z) {
				counter++;
			}
		}

		return -counter / 4.3363E-2;

	}

	private static double alphaHfderiv(SolutionR soln, int Z) {

		double sum = 0;

		for (int i = 0; i < soln.atoms.length; i++) {
			for (int j = i + 1; j < soln.atoms.length; j++) {
				sum += soln.atoms[i].crfParamDeriv(soln.atoms[j],
						getNum(soln.atomicNumbers[i], soln.atomicNumbers[j],
								Z));
			}
		}

		return sum / 4.3363E-2;
	}

	private static double betaderiv(int Z1, int Z2, int Z, int L1, int L2,
									int type) {

		double sum = 0;

		if (Z1 == Z && L1 == type) {
			sum += 0.5;
		}

		if (Z2 == Z && L2 == type) {
			sum += 0.5;
		}
		return sum;
	}

	private static int getNum(int Z1, int Z2, int Z) {
		int num = 0;

		if (Z1 == Z) {
			num += 1;
		}

		if (Z2 == Z) {
			num += 2;
		}

		return num - 1;
	}

	public static int index(NDDO6G orbital) {

		if (orbital.getL() == 0) {
			return 0;
		}
		else if (orbital.geti() == 1) {
			return 1;
		}
		else if (orbital.getj() == 1) {
			return 2;
		}
		else if (orbital.getk() == 1) {
			return 3;
		}

		return -1;
	}

	public static DoubleMatrix responseMatrix(SolutionR soln,
											  DoubleMatrix densityMatrixDeriv) {

		DoubleMatrix responsematrix =
				DoubleMatrix.zeros(soln.orbitals.length, soln.orbitals.length);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = j; k < soln.orbitals.length; k++) {
				double val = 0;
				if (j == k) {

					for (int l : soln.index[soln.atomNumber[j]]) {
						if (l > -1) {
							val += densityMatrixDeriv.get(l, l) *
									integralArray[integralcount];
							integralcount++;
						}
					}

					for (int l : soln.missingIndex[soln.atomNumber[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingIndex[soln.atomNumber[j]]) {
								if (m > -1) {
									if (soln.atomNumber[l] ==
											soln.atomNumber[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else if (soln.atomNumber[j] == soln.atomNumber[k]) {
					val += densityMatrixDeriv.get(j, k) *
							integralArray[integralcount];
					integralcount++;

					for (int l : soln.missingIndex[soln.atomNumber[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingIndex[soln.atomNumber[j]]) {
								if (m > -1) {
									if (soln.atomNumber[l] ==
											soln.atomNumber[m]) {
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
					for (int l : soln.index[soln.atomNumber[j]]) {
						if (l > -1) {
							for (int m : soln.index[soln.atomNumber[k]]) {
								if (m > -1) {
									val += densityMatrixDeriv.get(l, m) *
											integralArray[integralcount];
									integralcount++;
								}
							}
						}
					}
				}

				responsematrix.put(j, k, val);
				responsematrix.put(k, j, val);
			}
		}

		return responsematrix;
	}


	private static DoubleMatrix computeResponseVectorsLimited(DoubleMatrix x,
															  SolutionR soln) {//todo
		// duplicate from GeometrySecondDerivative

		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		DoubleMatrix densityMatrixDeriv =
				DoubleMatrix.zeros(soln.orbitals.length, soln.orbitals.length);

		for (int u = 0; u < densityMatrixDeriv.rows; u++) {
			for (int v = 0; v < densityMatrixDeriv.columns; v++) {
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

				densityMatrixDeriv.put(u, v, sum);
			}
		}

		DoubleMatrix responsematrix =
				DoubleMatrix.zeros(soln.orbitals.length, soln.orbitals.length);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = j; k < soln.orbitals.length; k++) {
				double val = 0;
				if (j == k) {

					for (int l : soln.index[soln.atomNumber[j]]) {
						if (l > -1) {
							val += densityMatrixDeriv.get(l, l) *
									integralArray[integralcount];
							integralcount++;
						}
					}

					for (int l : soln.missingIndex[soln.atomNumber[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingIndex[soln.atomNumber[j]]) {
								if (m > -1) {
									if (soln.atomNumber[l] ==
											soln.atomNumber[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else if (soln.atomNumber[j] == soln.atomNumber[k]) {
					val += densityMatrixDeriv.get(j, k) *
							integralArray[integralcount];
					integralcount++;

					for (int l : soln.missingIndex[soln.atomNumber[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingIndex[soln.atomNumber[j]]) {
								if (m > -1) {
									if (soln.atomNumber[l] ==
											soln.atomNumber[m]) {
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
					for (int l : soln.index[soln.atomNumber[j]]) {
						if (l > -1) {
							for (int m : soln.index[soln.atomNumber[k]]) {
								if (m > -1) {
									val += densityMatrixDeriv.get(l, m) *
											integralArray[integralcount];
									integralcount++;
								}
							}
						}
					}
				}

				responsematrix.put(j, k, val);
				responsematrix.put(k, j, val);
			}
		}

		DoubleMatrix R = DoubleMatrix.zeros(NOcc * NVirt, 1);

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


				R.put(count1, 0, element);

				count1++;
			}
		}

		DoubleMatrix p = new DoubleMatrix(NOcc * NVirt, 1);

		int counter = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				p.put(counter, 0, -R.get(counter, 0) +
						(soln.E.get(j + NOcc) - soln.E.get(i)) *
								x.get(counter));
				counter++;
			}
		}


		return p;
	}

	// fockDerivStatic is 1x8
	public static DoubleMatrix[] xArrayLimitedPople(SolutionR soln,
													DoubleMatrix[] fockDerivStaticPadded) {
		return xArrayLimitedThiel(soln, fockDerivStaticPadded);
//		int size = numNotNull(fockDerivStaticPadded);
//		DoubleMatrix[] fockDerivStatic = new DoubleMatrix[size];
//		size = 0;
//		for (DoubleMatrix dm : fockDerivStaticPadded) {
//			if (dm != null) {
//				fockDerivStatic[size] = dm.dup();
//				size++;
//			}
//		}
//		if (fockDerivStatic.length == 0)
//			return new DoubleMatrix[fockDerivStaticPadded.length];
//
//		int NOcc = (int) (soln.nElectrons / 2.0);
//		int NVirt = soln.orbitals.length - NOcc;
//
//		DoubleMatrix[] xArray = new DoubleMatrix[fockDerivStatic.length];
//		DoubleMatrix[] barray = new DoubleMatrix[fockDerivStatic.length];
//		DoubleMatrix[] parray = new DoubleMatrix[fockDerivStatic.length];
//		DoubleMatrix[] Farray = new DoubleMatrix[fockDerivStatic.length];
//		DoubleMatrix[] rarray = new DoubleMatrix[fockDerivStatic.length];
//
//		DoubleMatrix preconditioner = DoubleMatrix.zeros(NOcc * NVirt, 1);
//		DoubleMatrix preconditionerinv = DoubleMatrix.zeros(NOcc * NVirt, 1);
//
//		int counter = 0;
//
//		for (int i = 0; i < NOcc; i++) {
//			for (int j = 0; j < NVirt; j++) {
//
//				double e = (-soln.E.get(i) + soln.E.get(NOcc + j));
//
//				preconditioner.put(counter, Math.pow(e, -0.5));
//
//				preconditionerinv.put(counter, Math.pow(e, 0.5));
//
//				counter++;
//
//			}
//		}
//
//		DoubleMatrix D = DoubleMatrix.diag(preconditioner);
//
//		DoubleMatrix Dinv = DoubleMatrix.diag(preconditionerinv);
//
//		for (int a = 0; a < xArray.length; a++) {
//			DoubleMatrix F = DoubleMatrix.zeros(NOcc * NVirt, 1);
//
//			int count1 = 0;
//
//			for (int i = 0; i < NOcc; i++) {
//				for (int j = 0; j < NVirt; j++) {
//					double element = 0;
//
//					for (int u = 0; u < soln.orbitals.length; u++) {
//						for (int v = 0; v < soln.orbitals.length; v++) {
//							element +=
//									soln.C.get(i, u) * soln.C.get(j + NOcc,
//											v) *
//											fockDerivStatic[a].get(u, v);
//						}
//					}
//
//					element = element / (soln.E.get(j + NOcc) - soln.E.get(i));
//					F.put(count1, 0, element);
//					count1++;
//				}
//			}
//
//			F = D.mmul(F);
//			xArray[a] = DoubleMatrix.zeros(NOcc * NVirt, 1);
//			rarray[a] = xArray[a].dup();
//			barray[a] = F.dup();
//			Farray[a] = F.dup();
//		}
//
//
//		if (barray[0].rows == 0) {
//			DoubleMatrix[] densityderivs =
//					new DoubleMatrix[fockDerivStatic.length];
//
//			for (int i = 0; i < densityderivs.length; i++) {
//				densityderivs[i] = DoubleMatrix.zeros(0, 0);
//			}
//
//			return densityderivs;
//		}
//
//
//		ArrayList<DoubleMatrix> b = new ArrayList<>();
//
//		ArrayList<DoubleMatrix> p = new ArrayList<>();
//
//		int[] iterable = new int[barray.length];
//
//
//		DoubleMatrix F = DoubleMatrix.zeros(NOcc * NVirt, Farray.length);
//		for (int i = 0; i < Farray.length; i++) {
//			F.putColumn(i, Farray[i]);
//		}
//
//		while (GeometrySecondDerivative.numIterable(iterable) > 0) {
//
//			GeometrySecondDerivative.orthogonalise(barray);
//
////            System.err.println("only " + GeometrySecondDerivative
////            .numIterable(iterable) +
////            " left to go!");
//
//			for (int i = 0; i < barray.length; i++) {
//
//				b.add(barray[i].dup());
//				parray[i] = D.mmul(GeometrySecondDerivative
//						.computeResponseVectorsPople(Dinv.mmul(barray[i].dup()),
//								soln));
//				p.add(parray[i].dup());
//			}
//
//			for (int i = 0; i < barray.length; i++) {
//				DoubleMatrix newb = parray[i];
//
//				for (int j = 0; j < b.size(); j++) {
//					double num = b.get(j).transpose().mmul(parray[i]).get(0) /
//							b.get(j).transpose().mmul(b.get(j)).get(0);
//					newb = newb.sub(b.get(j).mmul(num));
//				}
//
//				barray[i] = newb.dup();
//			}
//
//			DoubleMatrix B = DoubleMatrix.zeros(NOcc * NVirt, b.size());
//			DoubleMatrix P = DoubleMatrix.zeros(NOcc * NVirt, b.size());
//
//			for (int i = 0; i < b.size(); i++) {
//				B.putColumn(i, b.get(i));
//
//				P.putColumn(i, b.get(i).sub(p.get(i)));
//			}
//
//
//			DoubleMatrix lhs = B.transpose().mmul(P);
//
//			DoubleMatrix rhs = B.transpose().mmul(F);
//			DoubleMatrix alpha;
//
//			try {
//				alpha = Solve.solve(lhs, rhs);
//			} catch (LapackException e) {
//				alpha = DoubleMatrix.ones(lhs.columns, rhs.columns);
////                System.err.println(Arrays.toString(fockDerivStatic));
////                System.err.println(lhs);
////                System.err.println(rhs);
//			}
//
//			for (int a = 0; a < xArray.length; a++) {
//
//				rarray[a] = DoubleMatrix.zeros(NOcc * NVirt, 1);
//				xArray[a] = DoubleMatrix.zeros(NOcc * NVirt, 1);
//			}
//
//			for (int i = 0; i < alpha.rows; i++) {
//				for (int j = 0; j < alpha.columns; j++) {
//
//					rarray[j] =
//							rarray[j].add((b.get(i).sub(p.get(i)))
//									.mmul(alpha.get(i,
//											j)));
//					xArray[j] = xArray[j].add(b.get(i).mmul(alpha.get(i, j)));
//				}
//			}
//
//			for (int j = 0; j < alpha.columns; j++) {
//				rarray[j] = rarray[j].sub(Farray[j]);
//
//				xArray[j] = Dinv.mmul(xArray[j]);
//
//				if (mag(rarray[j]) < 1E-10) { //todo play with this
//					iterable[j] = 1;
//				}
//				else if (Double.isNaN(mag(rarray[j]))) {
//					System.err.println(
//							"Pople algorithm fails; reverting to Thiel " +
//									"algorithm " +
//									"(don't" +
//									" " +
//									"panic)...");
//					return xArrayLimitedThiel(soln, fockDerivStaticPadded);
//				}
//				else {
//					iterable[j] = 0;
////                    System.err.println("convergence test: " + mag(rarray[j]));
//				}
//			}
//		}
//
//		DoubleMatrix[] xArrayPadded =
//				new DoubleMatrix[fockDerivStaticPadded.length];
//		size = 0;
//		for (int i = 0; i < fockDerivStaticPadded.length; i++) {
//			if (fockDerivStaticPadded[i] != null) {
//				xArrayPadded[i] = xArray[size];
//				size++;
//			}
//		}
//		return xArrayPadded;
	}

	private static DoubleMatrix[] xArrayLimitedThiel(SolutionR soln,
													 DoubleMatrix[] fockDerivStatic) {
		boolean bad = true;
		for (DoubleMatrix dm : fockDerivStatic) {
			if (dm != null) {
				bad = false;
				break;
			}
		}
		if (bad) return new DoubleMatrix[fockDerivStatic.length];
		int NOcc = (int) (soln.nElectrons / 2.0);
		int NVirt = soln.orbitals.length - NOcc;
		DoubleMatrix[] Farray = new DoubleMatrix[fockDerivStatic.length];
		DoubleMatrix[] xarray = new DoubleMatrix[fockDerivStatic.length];
		DoubleMatrix[] rarray = new DoubleMatrix[fockDerivStatic.length];
		DoubleMatrix[] dirs = new DoubleMatrix[fockDerivStatic.length];

		for (int a = 0; a < Farray.length; a++) {
			if (fockDerivStatic[a] != null) {
				DoubleMatrix F = DoubleMatrix.zeros(NOcc * NVirt, 1);

				int count1 = 0;

				for (int i = 0; i < NOcc; i++) {
					for (int j = 0; j < NVirt; j++) {

						double element = 0;

						for (int u = 0; u < soln.orbitals.length; u++) {
							for (int v = 0; v < soln.orbitals.length; v++) {
								element += soln.C.get(i, u) *
										soln.C.get(j + NOcc, v) *
										fockDerivStatic[a].get(u, v);
							}
						}

						F.put(count1, 0, element);

						count1++;
					}
				}

				Farray[a] = F;
				xarray[a] = DoubleMatrix.zeros(NOcc * NVirt, 1);
				rarray[a] = F.dup();
				dirs[a] = F.dup();
			}
		}

		// todo adrian what the frick is this
		for (DoubleMatrix dir : dirs) {
			if (dir != null) {
				if (dir.rows == 0) {
					DoubleMatrix[] densityderivs =
							new DoubleMatrix[fockDerivStatic.length];

					for (int i = 0; i < densityderivs.length; i++) {
						densityderivs[i] = DoubleMatrix.zeros(0, 0);
					}

					return densityderivs;
				}
			}
		}


		while (numNotNull(rarray) > 0) {
//            System.err.println("It's still running, don't worry: " +
//            numNotNull
			//            (rarray));
			ArrayList<DoubleMatrix> d = new ArrayList<>();
			for (int a = 0; a < rarray.length; a++) {
				if (rarray[a] != null) d.add(dirs[a].dup());
			}

			ArrayList<DoubleMatrix> p = new ArrayList<>();

			for (DoubleMatrix doubleMatrix : d) {
				p.add(computeResponseVectorsLimited(doubleMatrix, soln));
			}

			DoubleMatrix solver = new DoubleMatrix(p.size(), p.size());

			DoubleMatrix rhsvec = DoubleMatrix.zeros(p.size(), rarray.length);

			for (int a = 0; a < rhsvec.columns; a++) {
				if (rarray[a] != null) {
					DoubleMatrix rhs = new DoubleMatrix(p.size(), 1);

					for (int i = 0; i < rhs.rows; i++) {
						rhs.put(i, 0,
								2 * rarray[a].transpose().mmul(d.get(i)).get(0,
										0));

					}

					rhsvec.putColumn(a, rhs);
				}
			}

			for (int i = 0; i < solver.rows; i++) {
				for (int j = i; j < solver.rows; j++) {

					double val = p.get(j).transpose().mmul(d.get(i)).get(0,
							0) +
							p.get(i).transpose().mmul(d.get(j)).get(0, 0);
					solver.put(i, j, val);
					solver.put(j, i, val);
				}
			}

			DoubleMatrix alpha;

			try {
				alpha = Solve.solve(solver, rhsvec);
			} catch (LapackException e) {
				alpha = DoubleMatrix.ones(solver.columns, rhsvec.columns);
			}

			for (int a = 0; a < rhsvec.columns; a++) {
				if (rarray[a] != null) {


					for (int i = 0; i < alpha.rows; i++) {
						xarray[a] =
								xarray[a].add(d.get(i).mmul(alpha.get(i, a)));
						rarray[a] =
								rarray[a].sub(p.get(i).mmul(alpha.get(i, a)));

					}

					if (mag(rarray[a]) < 1E-7) {
						rarray[a] = null;
					}
				}
			}


			solver = new DoubleMatrix(solver.rows, solver.rows);

			for (int a = 0; a < rhsvec.columns; a++) {
				if (rarray[a] != null) {
					DoubleMatrix rhs = new DoubleMatrix(solver.rows, 1);

					for (int i = 0; i < rhs.rows; i++) {
						rhs.put(i, 0, -rarray[a].transpose().mmul(p.get(i))
								.get(0, 0));

					}

					rhsvec.putColumn(a, rhs);
				}
			}


			for (int i = 0; i < solver.rows; i++) {
				for (int j = 0; j < solver.rows; j++) {
					solver.put(i, j,
							d.get(j).transpose().mmul(p.get(i)).get(0, 0));
				}
			}

			DoubleMatrix beta;

			try {
				beta = Solve.solve(solver, rhsvec);
			} catch (Exception e) {
				beta = DoubleMatrix.ones(solver.columns, rhsvec.columns);
			}

			for (int a = 0; a < rhsvec.columns; a++) {
				if (rarray[a] != null) {


					dirs[a] = rarray[a].dup();

					for (int i = 0; i < beta.rows; i++) {
						dirs[a] = dirs[a].add(d.get(i).mmul(beta.get(i, a)));
					}
				}
			}


		}
		return xarray;
	}


	public static DoubleMatrix xArrayComplementary(SolutionR soln,
												   DoubleMatrix fockderiv) {
		int NOcc = (int) (soln.nElectrons / 2.0);
		if (NOcc == 0) {
			return DoubleMatrix.zeros(0, 0);
		}
		DoubleMatrix x = DoubleMatrix.zeros(NOcc - 1, 1);
		int count1 = 0;

		for (int j = 0; j < NOcc - 1; j++) {
			double element = 0;

			for (int u = 0; u < soln.orbitals.length; u++) {
				for (int v = 0; v < soln.orbitals.length; v++) {
					element += soln.C.get(NOcc - 1, u) * soln.C.get(j, v) *
							fockderiv.get(u, v);
				}
			}
			x.put(count1, 0,
					-element / (soln.E.get(NOcc - 1) - soln.E.get(j)));
			count1++;
		}
		return x;
	}


	public static DoubleMatrix densityDerivativeLimited(SolutionR soln,
														DoubleMatrix x) {

		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		DoubleMatrix densityMatrixDeriv =
				DoubleMatrix.zeros(soln.orbitals.length, soln.orbitals.length);
		for (int u = 0; u < densityMatrixDeriv.rows; u++) {
			for (int v = 0; v < densityMatrixDeriv.columns; v++) {
				double sum = 0;
				int count = 0;
				for (int i = 0; i < NOcc; i++) {
					for (int j = 0; j < NVirt; j++) {
							sum -= 2 * (soln.C.get(i, u) * soln.C.get(j + NOcc,
									v) +
									soln.C.get(j + NOcc, u) *
											soln.C.get(i, v)) *
									x.get(count, 0);
							count++;
						}
					}
				densityMatrixDeriv.put(u, v, sum);
			}
		}

		return densityMatrixDeriv;
	}

	public static DoubleMatrix xarrayForIE(SolutionR soln,
										   DoubleMatrix xlimited,
										   DoubleMatrix xcomplementary) {

		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		DoubleMatrix x = new DoubleMatrix(soln.orbitals.length - 1, 1);

		int count = 0;

		for (int i = xlimited.length - NVirt; i < xlimited.length; i++) {
			if (i > -1) {
				x.put(count, 0, xlimited.get(i, 0));
				count++;
			}
		}

		for (int i = 0; i < xcomplementary.length; i++) {
			x.put(count, 0, xcomplementary.get(i, 0));
			count++;
		}

		return x;
	}


	public static DoubleMatrix HOMOCoefficientDerivativeComplementary(
			DoubleMatrix x,
			SolutionR soln) {


		DoubleMatrix CDeriv = DoubleMatrix.zeros(1, soln.orbitals.length);

		int NOcc = (int) (soln.nElectrons / 2.0);

		for (int u = 0; u < soln.orbitals.length; u++) {
			double sum = 0;

			for (int k = 0; k < soln.orbitals.length; k++) {

				if (k < NOcc - 1) {
					sum -= soln.C.get(k, u) * x.get(k, 0);
				}
				else if (k >= NOcc) {
					if (k > 0) {
						sum -= soln.C.get(k, u) * x.get(k - 1, 0);
					}
				}
			}

			CDeriv.put(0, u, sum);
		}

		return CDeriv;
	}


	public static double MNDOIEDeriv(SolutionR soln,
									 DoubleMatrix coeffDeriv,
									 DoubleMatrix Fderiv) {

		int index = (int) (soln.nElectrons / 2.0) - 1;

		if (index < 0) {
			return 0;
		}

		DoubleMatrix coeff = soln.C.getRow(index);

		double sum = 0;

		for (int i = 0; i < soln.orbitals.length; i++) {
			for (int j = 0; j < soln.orbitals.length; j++) {
				sum += soln.F.get(i, j) * (coeff.get(i) * coeffDeriv.get(j) +
						coeff.get(j) * coeffDeriv.get(i));
				sum += Fderiv.get(i, j) * coeff.get(i) * coeff.get(j);
			}
		}

		return sum;


	}

	public static double MNDODipoleDeriv(SolutionR soln,
										 DoubleMatrix densityderiv, int Z,
										 int paramnum) {

		if (soln.dipole <= 1E-6) {
			return 0;
		}

		NDDOAtom[] atoms = soln.atoms;

		double D1deriv = 0;

		if (paramnum == 5 || paramnum == 6) {
			for (int i = 0; i < atoms.length; i++) {
				if (atoms[i].getAtomProperties().getZ() == Z) {
					D1deriv = D1Deriv(atoms[i], paramnum - 5);
					break;
				}
			}
		}

		int[][] index = soln.index;

		DoubleMatrix densityMatrix = soln.densityMatrix();

		double[] populationsderiv = new double[atoms.length];

		for (int j = 0; j < atoms.length; j++) {
			double sum = 0;
			for (int k : index[j]) {
				if (k > -1) {
					sum += densityderiv.get(k, k);
				}
			}

			populationsderiv[j] = -sum;
		}


		double[] com = new double[]{0, 0, 0};

		double mass = 0;

		for (NDDOAtom atom : atoms) {
			com[0] = com[0] + atom.getMass() * atom.getCoordinates()[0];
			com[1] = com[1] + atom.getMass() * atom.getCoordinates()[1];
			com[2] = com[2] + atom.getMass() * atom.getCoordinates()[2];
			mass += atom.getMass();
		}

		com[0] = com[0] / mass;
		com[1] = com[1] / mass;
		com[2] = com[2] / mass;


		double[] chargedip = new double[]{0, 0, 0};

		for (int j = 0; j < atoms.length; j++) {
			chargedip[0] += 2.5416 * populationsderiv[j] *
					(atoms[j].getCoordinates()[0] - com[0]);
			chargedip[1] += 2.5416 * populationsderiv[j] *
					(atoms[j].getCoordinates()[1] - com[1]);
			chargedip[2] += 2.5416 * populationsderiv[j] *
					(atoms[j].getCoordinates()[2] - com[2]);
		}


		double[] hybridip = new double[]{0, 0, 0};

		for (int j = 0; j < atoms.length; j++) {

			if (index[j][1] != -1) {//exclude hydrogen
				hybridip[0] = hybridip[0] - 2.5416 * 2 * atoms[j].D1 *
						densityderiv.get(index[j][0], index[j][1]);
				hybridip[1] = hybridip[1] - 2.5416 * 2 * atoms[j].D1 *
						densityderiv.get(index[j][0], index[j][2]);
				hybridip[2] = hybridip[2] - 2.5416 * 2 * atoms[j].D1 *
						densityderiv.get(index[j][0], index[j][3]);

				if (atoms[j].getAtomProperties().getZ() == Z) {
					hybridip[0] = hybridip[0] - 2.5416 * 2 * D1deriv *
							densityMatrix.get(index[j][0], index[j][1]);
					hybridip[1] = hybridip[1] - 2.5416 * 2 * D1deriv *
							densityMatrix.get(index[j][0], index[j][2]);
					hybridip[2] = hybridip[2] - 2.5416 * 2 * D1deriv *
							densityMatrix.get(index[j][0], index[j][3]);
				}
			}


		}


		double[] dipoletot =
				new double[]{chargedip[0] + hybridip[0],
						chargedip[1] + hybridip[1],
						chargedip[2] + hybridip[2]};


		return (dipoletot[0] * soln.dipoletot[0] +
				dipoletot[1] * soln.dipoletot[1] +
				dipoletot[2] * soln.dipoletot[2]) / soln.dipole;
	}


	private static double mag(DoubleMatrix gradient) {

		double sum = 0;
		for (int i = 0; i < gradient.length; i++) {
			sum += gradient.get(i) * gradient.get(i);
		}

		return Math.sqrt(sum);
	}

	private static int numNotNull(DoubleMatrix[] rarray) {

		int count = 0;
		for (DoubleMatrix r : rarray) {
			if (r != null) {
				count++;
			}
		}

		return count;
	}


	private static double D1Derivfinite(NDDOAtom a, int type) throws Exception {

		double D1 = a.D1;


		NDDOParams params = a.getParams().clone();
		params.modifyParam(5 + type, Utils.LAMBDA);

		Class<? extends NDDOAtom> c = a.getClass();
		Constructor ctor =
				c.getDeclaredConstructor(c, a.getParams().getClass());
		ctor.setAccessible(true);

		NDDOAtom a2 = (NDDOAtom) ctor.newInstance(a, params);

		double D1perturbed = a2.D1;

		return (D1perturbed - D1) / Utils.LAMBDA;

	}

	private static double D2Derivfinite(NDDOAtom a, int type) throws Exception {

		double D2 = a.D2;


		NDDOParams params = a.getParams().clone();
		params.modifyParam(5 + type, Utils.LAMBDA);

		Class<? extends NDDOAtom> c = a.getClass();
		Constructor ctor =
				c.getDeclaredConstructor(c, a.getParams().getClass());
		ctor.setAccessible(true);

		NDDOAtom a2 = (NDDOAtom) ctor.newInstance(a, params);

		double D2perturbed = a2.D2;

		return (D2perturbed - D2) / Utils.LAMBDA;

	}

	private static double p1Derivfinite(NDDOAtom a, int type) throws Exception {

		double p1 = a.p1;

		NDDOParams params = a.getParams().clone();
		params.modifyParam(5 + type, Utils.LAMBDA);

		Class<? extends NDDOAtom> c = a.getClass();
		Constructor ctor =
				c.getDeclaredConstructor(c, a.getParams().getClass());
		ctor.setAccessible(true);

		NDDOAtom a2 = (NDDOAtom) ctor.newInstance(a, params);

		double p1perturbed = a2.p1;

		return (p1perturbed - p1) / Utils.LAMBDA;

	}

	private static double p2Derivfinite(NDDOAtom a, int type) throws Exception {

		double p2 = a.p2;

		NDDOParams params = a.getParams().clone();
		params.modifyParam(5 + type, Utils.LAMBDA);

		Class<? extends NDDOAtom> c = a.getClass();
		Constructor ctor =
				c.getDeclaredConstructor(c, a.getParams().getClass());
		ctor.setAccessible(true);

		NDDOAtom a2 = (NDDOAtom) ctor.newInstance(a, params);

		double p2perturbed = a2.p2;

		return (p2perturbed - p2) / Utils.LAMBDA;

	}


}
