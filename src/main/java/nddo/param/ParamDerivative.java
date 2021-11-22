package nddo.param;

import nddo.NDDO6G;
import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.solution.SolutionR;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;
import scf.GTO;
import scf.LCGTO;
import tools.Utils;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Arrays;

import static nddo.geometry.GeometrySecondDerivative.*;

public class ParamDerivative {

	private static double p1Deriv(NDDOAtom a, int type) {

		double D1deriv = D1Deriv(a, type);

		if (a.getAtomProperties().getZ() == 1) {
			return 0;
		}

		return -a.p1 * a.p1 * a.D1 / (a.p1 * a.p1 * a.p1 -
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
							 double D21, double p02, double p12, double p22,
							 double D12, double D22, double R, int num,
							 double D1deriv, double D2deriv, double p1deriv,
							 double p2deriv) {
		return 0;
	}

	private static double quz(double p01, double p11, double p21, double D11,
							  double D21, double p02, double p12, double p22,
							  double D12, double D22, double R, int num,
							  double D1deriv, double D2deriv, double p1deriv,
							  double p2deriv) {
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
							   double D21, double p02, double p12, double p22,
							   double D12, double D22, double R, int num,
							   double D1deriv, double D2deriv, double p1deriv,
							   double p2deriv) {
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
							   double D21, double p02, double p12, double p22,
							   double D12, double D22, double R, int num,
							   double D1deriv, double D2deriv, double p1deriv,
							   double p2deriv) {
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
								  double D11, double D21, double p02,
								  double p12, double p22, double D12,
								  double D22, double R, int num,
								  double D1deriv,
								  double D2deriv, double p1deriv,
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
					Math.pow((R - D22) * (R - D22) + (D11 - D22) * (D11 - D22) +
							a12 * a12, -1.5)
					- 0.25 * ((2 * D22 + D11 - R) * D2deriv + a12 * p2deriv) *
					Math.pow((R - D22) * (R - D22) + (D11 + D22) * (D11 + D22) +
							a12 * a12, -1.5)
					- 0.25 * ((2 * D22 + R - D11) * D2deriv + a12 * p2deriv) *
					Math.pow((R + D22) * (R + D22) + (D11 - D22) * (D11 - D22) +
							a12 * a12, -1.5)
					+ 0.25 * ((2 * D22 + D11 + R) * D2deriv + a12 * p2deriv) *
					Math.pow((R + D22) * (R + D22) + (D11 + D22) * (D11 + D22) +
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
								  double D11, double D21, double p02,
								  double p12, double p22, double D12,
								  double D22, double R, int num,
								  double D1deriv,
								  double D2deriv, double p1deriv,
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
								double D1deriv, double D2deriv, double p1deriv,
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
					a12 * (p1deriv + p2deriv)) * Math.pow(
					(R + D11 - 2 * D21) * (R + D11 - 2 * D21) + a12 * a12,
					-1.5)
					- 0.125 * ((D11 + 2 * D21 - R) * (D1deriv + 2 * D2deriv) +
					a12 * (p1deriv + p2deriv)) * Math.pow(
					(R - D11 - 2 * D21) * (R - D11 - 2 * D21) + a12 * a12,
					-1.5)
					+ 0.125 * ((R + D11 + 2 * D21) * (D1deriv + 2 * D2deriv) +
					a12 * (p1deriv + p2deriv)) * Math.pow(
					(R + D11 + 2 * D21) * (R + D11 + 2 * D21) + a12 * a12,
					-1.5)
					- 0.125 * ((R - D11 + 2 * D21) * (2 * D2deriv - D1deriv) +
					a12 * (p1deriv + p2deriv)) * Math.pow(
					(R - D11 + 2 * D21) * (R - D11 + 2 * D21) + a12 * a12,
					-1.5)
					- 0.25 * ((R + D11) * D1deriv + a12 * (p1deriv + p2deriv)) *
					Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
					+ 0.25 * ((D11 - R) * D1deriv + a12 * (p1deriv + p2deriv)) *
					Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
		}

		return 0;

	}

	private static double QpipiQpipi(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double R, int num,
									 double D1deriv, double D2deriv,
									 double p1deriv, double p2deriv) {
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
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R, int num,
								   double D1deriv, double D2deriv,
								   double p1deriv, double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) {
			return -0.125 * (4 * D21 * D2deriv + a22 * p2deriv) * Math.pow(
					(R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 + a22 * a22,
					-1.5)
					- 0.125 * (4 * D21 * D2deriv + a22 * p2deriv) * Math.pow(
					(R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 + a22 * a22,
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
							a22 * a22, -1.5)
					- 0.125 * (2 * (2 * D22 + R) * D2deriv + a22 * p2deriv) *
					Math.pow((R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 +
							a22 * a22, -1.5)
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
			return -0.125 * ((2 * (2 * D21 - R) + 4 * D21) * D2deriv +
					2 * a22 * p2deriv) * Math.pow(
					(R - 2 * D21) * (R - 2 * D21) + 4 * D21 * D21 + a22 * a22,
					-1.5)
					- 0.125 * ((2 * (2 * D21 + R) + 4 * D21) * D2deriv +
					2 * a22 * p2deriv) * Math.pow(
					(R + 2 * D21) * (R + 2 * D21) + 4 * D21 * D21 + a22 * a22,
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
					Math.pow((R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) +
							a22 * a22, -1.5)
					- 0.0625 *
					(2 * (R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) +
							a22 * a22, -1.5)
					- 0.0625 *
					(2 * (2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) +
							a22 * a22, -1.5)
					- 0.0625 *
					(2 * (2 * D21 - R - 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) +
							a22 * a22, -1.5)
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
					Math.pow((R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) +
							a22 * a22, -1.5)
					- 0.0625 *
					(2 * (R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) +
							a22 * a22, -1.5)
					- 0.0625 *
					(2 * (2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) +
							a22 * a22, -1.5)
					- 0.0625 *
					(2 * (R - 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) +
							a22 * a22, -1.5)
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
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R, int num,
								   double D1deriv, double D2deriv,
								   double p1deriv, double p2deriv) {
		double a22 = p21 + p22;
		if (num == 0) {
			return -0.125 *
					((R + 2 * D21 - 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 - D22) * (R + D21 - D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					+ 0.125 * ((R + 2 * D21) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 - D22) * (R + D21 - D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					+ 0.125 * ((R + 2 * D21) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 + D22) * (R + D21 + D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					- 0.125 *
					((R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 + D22) * (R + D21 + D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					+ 0.125 * ((2 * D21 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 - D22) * (R - D21 - D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					- 0.125 *
					((2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 - D22) * (R - D21 - D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					- 0.125 *
					((2 * D21 - 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 + D22) * (R - D21 + D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					+ 0.125 * ((2 * D21 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 + D22) * (R - D21 + D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5);
		}
		if (num == 1) {
			return -0.125 *
					((2 * D22 - 2 * D21 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 - D22) * (R + D21 - D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					+ 0.125 * ((2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 - D22) * (R + D21 - D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					+ 0.125 * ((R + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 + D22) * (R + D21 + D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					- 0.125 *
					((R + 2 * D21 + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R + D21 + D22) * (R + D21 + D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					+ 0.125 * ((2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 - D22) * (R - D21 - D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					- 0.125 *
					((2 * D21 + 2 * D22 - R) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 - D22) * (R - D21 - D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
					- 0.125 *
					((2 * D22 + R - 2 * D21) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 + D22) * (R - D21 + D22) +
							(D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
					+ 0.125 * ((R + 2 * D22) * D2deriv + a22 * p2deriv) *
					Math.pow((R - D21 + D22) * (R - D21 + D22) +
							(D21 + D22) * (D21 + D22) + a22 * a22, -1.5);
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
					- 0.125 * ((2 * (2 * D21 + R) + 4 * D21) * D2deriv +
					2 * a22 * p2deriv) * Math.pow(
					(R + 2 * D21) * (R + 2 * D21) + 4 * D21 * D21 + a22 * a22,
					-1.5)
					- 0.125 * ((2 * (2 * D21 - R) + 4 * D21) * D2deriv +
					2 * a22 * p2deriv) * Math.pow(
					(R - 2 * D21) * (R - 2 * D21) + 4 * D21 * D21 + a22 * a22,
					-1.5);


		}

		return 0;
	}

	private static double ssssderiv(double p01, double p11, double p21,
									double D11, double D21, double p02,
									double p12, double p22, double D12,
									double D22, double R, int num,
									double D1deriv, double D2deriv,
									double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double ssppippideriv(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double R, int num,
										double D1deriv, double D2deriv,
										double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double sspzpzderiv(double p01, double p11, double p21,
									  double D11, double D21, double p02,
									  double p12, double p22, double D12,
									  double D22, double R, int num,
									  double D1deriv, double D2deriv,
									  double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
						D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double ppippissderiv(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double R, int num,
										double D1deriv, double D2deriv,
										double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pzpzssderiv(double p01, double p11, double p21,
									  double D11, double D21, double p02,
									  double p12, double p22, double D12,
									  double D22, double R, int num,
									  double D1deriv, double D2deriv,
									  double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double ppippippippideriv(double p01, double p11, double p21,
											double D11, double D21, double p02,
											double p12, double p22, double D12,
											double D22, double R, int num,
											double D1deriv, double D2deriv,
											double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				QpipiQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pxpxpypyderiv(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double R, int num,
										double D1deriv, double D2deriv,
										double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				QxxQyy(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double ppippipzpzderiv(double p01, double p11, double p21,
										  double D11, double D21, double p02,
										  double p12, double p22, double D12,
										  double D22, double R, int num,
										  double D1deriv, double D2deriv,
										  double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
						D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				QpipiQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pzpzppippideriv(double p01, double p11, double p21,
										  double D11, double D21, double p02,
										  double p12, double p22, double D12,
										  double D22, double R, int num,
										  double D1deriv, double D2deriv,
										  double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				QpipiQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pzpzpzpzderiv(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double R, int num,
										double D1deriv, double D2deriv,
										double p1deriv, double p2deriv) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
						D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				QzzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double spzssderiv(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double R, int num,
									 double D1deriv, double D2deriv,
									 double p1deriv, double p2deriv) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
				f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double spzppippideriv(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double R, int num,
										 double D1deriv, double D2deriv,
										 double p1deriv, double p2deriv) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
				f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				uzQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
						num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double spzpzpzderiv(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double R, int num,
									   double D1deriv, double D2deriv,
									   double p1deriv, double p2deriv) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
				f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				uzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
						D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double ssspzderiv(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double R, int num,
									 double D1deriv, double D2deriv,
									 double p1deriv, double p2deriv) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double ppippispzderiv(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double R, int num,
										 double D1deriv, double D2deriv,
										 double p1deriv, double p2deriv) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) -
				uzQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pzpzspzderiv(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double R, int num,
									   double D1deriv, double D2deriv,
									   double p1deriv, double p2deriv) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv) -
				uzQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
						f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double sppisppideriv(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double R, int num,
										double D1deriv, double D2deriv,
										double p1deriv, double p2deriv) {
		return upiupi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double spzspzderiv(double p01, double p11, double p21,
									  double D11, double D21, double p02,
									  double p12, double p22, double D12,
									  double D22, double R, int num,
									  double D1deriv, double D2deriv,
									  double p1deriv, double p2deriv) {
		return uzuz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num,
				D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double sppippipzderiv(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double R, int num,
										 double D1deriv, double D2deriv,
										 double p1deriv, double p2deriv) {
		return upiQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
				num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double ppipzsppideriv(double p01, double p11, double p21,
										 double D11, double D21, double p02,
										 double p12, double p22, double D12,
										 double D22, double R, int num,
										 double D1deriv, double D2deriv,
										 double p1deriv, double p2deriv) {
		return -upiQpiz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R,
				f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double ppipzppipzderiv(double p01, double p11, double p21,
										  double D11, double D21, double p02,
										  double p12, double p22, double D12,
										  double D22, double R, int num,
										  double D1deriv, double D2deriv,
										  double p1deriv, double p2deriv) {
		return QpizQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R,
				num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	private static double pxpypxpyderiv(double p01, double p11, double p21,
										double D11, double D21, double p02,
										double p12, double p22, double D12,
										double D22, double R, int num,
										double D1deriv, double D2deriv,
										double p1deriv, double p2deriv) {
		return 0.5 *
				(ppippippippideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12,
						D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) -
						pxpxpypyderiv(p01, p11, p21, D11, D21, p02, p12, p22,
								D12, D22, R, num, D1deriv, D2deriv, p1deriv,
								p2deriv));
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
												   NDDO6G d, double D1deriv,
												   double D2deriv,
												   double p1deriv,
												   double p2deriv, int num,
												   int type) {


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
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2
												, R, num, D1deriv, D2deriv,
												p1deriv, p2deriv);

									case 1:
										if (d.getk() == 1) {//(ss|spz)
											return ssspzderiv(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
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
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);

										case 1:
											if (d.getk() == 1) {//(ss|pzpz)
												return sspzpzderiv(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);
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
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2
												, R, num, D1deriv, D2deriv,
												p1deriv, p2deriv);
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
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);

										case 1:
											if (d.getk() == 1) {//(spz|spz)
												return spzspzderiv(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);
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
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);

											case 1:
												if (d.getk() == 1) {//(spz
													// |pzpz)
													return spzpzpzderiv(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
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
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
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
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2
												, R, num, D1deriv, D2deriv,
												p1deriv, p2deriv);
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
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
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
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
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
																, R, num,
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
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);

										case 1:
											if (d.getk() == 1) {//(pzs|spz)
												return spzspzderiv(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);
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
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);

											case 1:
												if (d.getk() == 1) {//(pzs
													// |pzpz)
													return spzpzpzderiv(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
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
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
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
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);

											case 1:
												if (d.getk() == 1) {//(pzpz
													// |spz)
													return pzpzspzderiv(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
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
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
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
																, R, num,
																D1deriv,
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
														a.p1, a.p2, a.D1, a.D2,
														c.p0, c.p1, c.p2, c.D1,
														c.D2, R, num, D1deriv,
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
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
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
														a.p1, a.p2, a.D1, a.D2,
														c.p0, c.p1, c.p2, c.D1,
														c.D2, R, num, D1deriv,
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
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzderiv(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	R, num,
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
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2
												, R, num, D1deriv, D2deriv,
												p1deriv, p2deriv);
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
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
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
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
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
																, R, num,
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
						case 1:
							if (b.getk() == 1) {//(ppipz|??)
								switch (c.getL()) {
									case 0://(ppipz|s?)
										if (d.geti() == a.geti() &&
												d.getj() == a.getj() &&
												d.getk() == 0) {//(ppipz|sppi)
											return ppipzsppideriv(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
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
														a.p1, a.p2, a.D1, a.D2,
														c.p0, c.p1, c.p2, c.D1,
														c.D2, R, num, D1deriv,
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
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzderiv(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	R, num,
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
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
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
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
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
																, R, num,
																D1deriv,
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
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R, num,
																D1deriv,
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
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, R,
																num, D1deriv,
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
																, R, num,
																D1deriv,
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
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
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
								   int num, int type) {
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
											D[l], D1deriv, D2deriv, p1deriv,
											p2deriv, num, type) * 27.21;
						}


					}
				}
			}
		}


		return sum2;
	}

@Deprecated
	public static double getGderivfinite(NDDO6G a, NDDO6G b, NDDO6G c,
										 NDDO6G d,
										 int num, int type) {

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
//				System.exit(0);
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
//				System.exit(0);
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
//				System.exit(0);
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
//				System.exit(0);
			}
		}

		double finalval =
				NDDO6G.getS(A.getOrbitals()[aindex], B.getOrbitals()[bindex]);

		return (finalval - initial) / Utils.LAMBDA;


	}

	@Deprecated
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
//				System.exit(0);
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
//				System.exit(0);
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
//		System.exit(0);
		return 0;

	}

	@Deprecated
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


	public static SimpleMatrix[][] MNDOStaticMatrixDeriv(SolutionR soln, int Z,
														 int firstParamIndex) {
		NDDOAtom[] atoms = soln.atoms;
		SimpleMatrix[] HDerivs = new SimpleMatrix[8];
		SimpleMatrix[] FDerivs = new SimpleMatrix[8];

		if (firstParamIndex <= 1) HDerivs[1] = betafockderivstatic(soln, Z, 0);
		if (firstParamIndex <= 1) FDerivs[1] = HDerivs[1].copy();
		if (firstParamIndex <= 3) HDerivs[3] = uxxfockderivstatic(soln, Z, 0);
		if (firstParamIndex <= 3) FDerivs[3] = HDerivs[3].copy();
		if (firstParamIndex <= 5)
			HDerivs[5] = zetaHderivstatic(atoms, soln, Z, 0);
		if (firstParamIndex <= 5)
			FDerivs[5] =
					HDerivs[5].copy().plus(zetaGderivstatic(atoms, soln, Z,
							0));

		if (Z != 1) {
			if (firstParamIndex <= 2)
				HDerivs[2] = betafockderivstatic(soln, Z, 1);
			if (firstParamIndex <= 2) FDerivs[2] = HDerivs[2].copy();
			if (firstParamIndex <= 4)
				HDerivs[4] = uxxfockderivstatic(soln, Z, 1);
			if (firstParamIndex <= 4) FDerivs[4] = HDerivs[4].copy();
			if (firstParamIndex <= 6)
				HDerivs[6] = zetaHderivstatic(atoms, soln, Z, 1);
			if (firstParamIndex <= 6)
				FDerivs[6] = HDerivs[6].copy()
						.plus(zetaGderivstatic(atoms, soln, Z, 1));
		}
		return new SimpleMatrix[][]{HDerivs, FDerivs};
	}

	public static double MNDOHFDeriv(SolutionR soln, SimpleMatrix Hderiv,
									 SimpleMatrix Fderiv) {

		double e = 0;

		SimpleMatrix densitymatrix = soln.densityMatrix();

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = 0; k < soln.orbitals.length; k++) {
				e += 0.5 * densitymatrix.get(j, k) *
						(Hderiv.get(j, k) + Fderiv.get(j, k));
			}
		}

		return e / 4.3363E-2;
	}

	private static double zetaHfderiv(SolutionR soln, int Z, int type) {

		SimpleMatrix densitymatrix = soln.densityMatrix();


		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.orbsOfAtom;

		int[][] missingIndex = soln.missingOfAtom;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);

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
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { // case 3
					double Huk = NDDO6G.betaparamderiv(orbitals[j],
							orbitals[k],
							getNum(atomicnumbers[atomNumber[j]],
									atomicnumbers[atomNumber[k]], Z), type);
					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}

		SimpleMatrix G = new SimpleMatrix(orbitals.length, orbitals.length);

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

				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		}

		SimpleMatrix F = H.copy().plus(G);

		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * densitymatrix.get(j, k) *
						(H.get(j, k) + F.get(j, k));
			}
		}

		return e / 4.3363E-2;

	}

	@Deprecated
	private static SimpleMatrix zetafockderivstatic(NDDOAtom[] atoms,
													SolutionR soln, int Z,
													int type) {


		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.orbsOfAtom;

		int[][] missingIndex = soln.missingOfAtom;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);

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
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { // case 3
					double Huk = NDDO6G.betaparamderiv(orbitals[j],
							orbitals[k],
							getNum(atomicnumbers[atomNumber[j]],
									atomicnumbers[atomNumber[k]], Z), type);
					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}

		SimpleMatrix G = new SimpleMatrix(orbitals.length, orbitals.length);

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

				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		}

		SimpleMatrix F = H.copy().plus(G);

		return F;
	}

	private static SimpleMatrix zetaHderivstatic(NDDOAtom[] atoms,
												 SolutionR soln, int Z,
												 int type) {


		NDDO6G[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);

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
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { // case 3
					double Huk = NDDO6G.betaparamderiv(orbitals[j],
							orbitals[k],
							getNum(atomicnumbers[atomNumber[j]],
									atomicnumbers[atomNumber[k]], Z), type);
					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}

		return H;
	}

	private static SimpleMatrix zetaGderivstatic(NDDOAtom[] atoms,
												 SolutionR soln, int Z,
												 int type) {


		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.orbsOfAtom;

		int[][] missingIndex = soln.missingOfAtom;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		SimpleMatrix G = new SimpleMatrix(orbitals.length, orbitals.length);

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

				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		}


		return G;
	}


	private static double uxxHfderiv(SolutionR soln, int Z, int type) {

		SimpleMatrix densitymatrix = soln.densityMatrix();


		NDDO6G[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

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

	private static SimpleMatrix uxxfockderivstatic(SolutionR soln, int Z,
												   int type) {

		SimpleMatrix F =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		for (int j = 0; j < soln.orbitals.length; j++) {
			if (soln.atomicNumbers[soln.atomOfOrb[j]] == Z &&
					soln.orbitals[j].getL() == type) {
				F.set(j, j, 1);
			}
		}

		return F;

	}

	private static double betaHfderiv(SolutionR soln, int Z, int type) {

		SimpleMatrix densitymatrix = soln.densityMatrix();


		NDDO6G[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

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

	private static SimpleMatrix betafockderivstatic(SolutionR soln, int Z,
													int type) {

		SimpleMatrix F =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		NDDO6G[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

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
					F.set(j, k, H);
					F.set(k, j, H);
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

	private static int index(NDDO6G orbital) {

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

	public static SimpleMatrix responseMatrix(SolutionR soln,
											  SimpleMatrix densityMatrixDeriv) {

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

		return responsematrix;
	}


	private static SimpleMatrix computeResponseVectorsLimited(SimpleMatrix x,
															  SolutionR soln) {//todo
		// duplicate from GeometrySecondDerivative

		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix densityMatrixDeriv =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
			for (int v = 0; v < densityMatrixDeriv.numCols(); v++) {
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

	public static SimpleMatrix[] xArrayLimitedPople(SolutionR soln,
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
		SimpleMatrix[] xarrayHold = new SimpleMatrix[fockderivstatic.length];
		SimpleMatrix[] barray = new SimpleMatrix[length];
		SimpleMatrix[] parray = new SimpleMatrix[length];
		SimpleMatrix[] Farray = new SimpleMatrix[length];
		SimpleMatrix[] rarray = new SimpleMatrix[length];

		// configure preconditioners
		SimpleMatrix D = new SimpleMatrix(nonv, nonv, DMatrixSparseCSC.class);
		SimpleMatrix Dinv =
				new SimpleMatrix(nonv, nonv, DMatrixSparseCSC.class);

		int counter = 0;
		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = (-soln.E.get(i) + soln.E.get(NOcc + j));

				D.set(counter, counter, Math.pow(e, -0.5));
				Dinv.set(counter, counter, Math.pow(e, 0.5));
				counter++;
			}
		}

		// convert AO to MO basis
		for (int a = 0; a < length; a++) {
			SimpleMatrix F = new SimpleMatrix(nonv, 1);

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

					element = element / (soln.E.get(j + NOcc) - soln.E.get(i));

					F.set(count, 0, element);

					count++;
				}
			}

			F = D.mult(F);

			barray[a] = F;
			Farray[a] = F;
		}

		// main loop
		int[] iterable = new int[length];
		ArrayList<SimpleMatrix> prevBs = new ArrayList<>();
		ArrayList<SimpleMatrix> prevPs = new ArrayList<>();

		// convert Farray into matrix form
		SimpleMatrix F = new SimpleMatrix(nonv, length);
		for (int i = 0; i < Farray.length; i++) {
			F.setColumn(i, 0, Farray[i].getDDRM().data);
		}

		double[] oldrMags = new double[rarray.length];
		Arrays.fill(oldrMags, 1);

		bigLoop:
		while (Utils.numIterable(iterable) > 0) {
			Utils.orthogonalise(barray);

			for (int i = 0; i < length; i++) {
				prevBs.add(barray[i]);

				parray[i] = D.mult(computeResponseVectorsPople(
						Dinv.mult(barray[i]), soln));

				prevPs.add(parray[i]);
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix newb = parray[i];

				// orthogonalize against all previous Bs
				for (SimpleMatrix prevB : prevBs) {
					double num = prevB.transpose().mult(parray[i]).get(0) /
							prevB.transpose().mult(prevB).get(0);

					newb = newb.minus(prevB.scale(num));
				}

				barray[i] = newb;
			}

			// convert prevBs and prevPs into matrix form, transposed
			SimpleMatrix Bt = new SimpleMatrix(prevBs.size(), nonv);
			SimpleMatrix P = new SimpleMatrix(nonv, prevPs.size());
			for (int i = 0; i < prevBs.size(); i++) {
				Bt.setRow(i, 0, prevBs.get(i).getDDRM().data);
				P.setColumn(i, 0,
						prevBs.get(i).minus(prevPs.get(i)).getDDRM().data);
			}

			SimpleMatrix lhs = Bt.mult(P);
			SimpleMatrix rhs = Bt.mult(F);
			SimpleMatrix alpha;
			try {
				alpha = lhs.solve(rhs);
			} catch (SingularMatrixException e) {
				alpha = Utils.filled(lhs.numCols(), rhs.numCols(), 1);
			}

			// reset r and x array
			for (int a = 0; a < length; a++) {
				rarray[a] = new SimpleMatrix(nonv, 1);
				xarray[a] = new SimpleMatrix(nonv, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {
					// B with tilde
					rarray[j] = rarray[j].plus(
							prevBs.get(i)
									.minus(prevPs.get(i))
									.scale(alpha.get(i, j))
					);

					xarray[j] = xarray[j].plus(
							prevBs.get(i)
									.scale(alpha.get(i, j))
					);
				}
			}

			int xarrayHoldNN = Utils.numNotNull(xarrayHold);
			for (int j = 0; j < alpha.numCols(); j++) {
				rarray[j] = rarray[j].minus(Farray[j]);
				xarray[j] = Dinv.mult(xarray[j]);

				double mag = Utils.mag(rarray[j]);
				if (mag > oldrMags[j] || mag != mag) {
					if (xarrayHoldNN == xarrayHold.length) {
						soln.getRm().getLogger().warn(
								"Slight numerical instability detected; " +
										"returning lower precision values. " +
										"rarray mag = {}",
								mag);
						xarray = xarrayHold;
						break bigLoop;
					}
					else {
						if (mag > oldrMags[j]) {
							soln.getRm().getLogger().error(
									"Numerical instability detected; " +
											"reverting to Thiel algorithm.");
							return xArrayLimitedThiel(soln, fockderivstatic);
						}
						if (mag != mag) {
							soln.getRm().getLogger()
									.error("Pople algorithm fails; " +
											"reverting to Thiel algorithm...");
							return xArrayLimitedThiel(soln, fockderivstatic);
						}
					}
				}
				else {
					if (mag < 1E-7) {
						xarrayHold[j] = xarray[j];
						if (mag < 1E-10) {
							iterable[j] = 1;
						}
					}
					else {
						iterable[j] = 0;
						soln.getRm().getLogger().trace(
								"Pople convergence test: " + mag);
					}
				}

				oldrMags[j] = mag;
			}
		}
		return xarray;
	}

	private static SimpleMatrix[] xArrayLimitedThiel(SolutionR soln,
													 SimpleMatrix[] fockDerivStatic) {
		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix[] xarray = new SimpleMatrix[fockDerivStatic.length];
		SimpleMatrix[] rarray = new SimpleMatrix[fockDerivStatic.length];
		SimpleMatrix[] dirs = new SimpleMatrix[fockDerivStatic.length];

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
									v) * fockDerivStatic[a].get(u, v);
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
					new SimpleMatrix[fockDerivStatic.length];

			for (int i = 0; i < densityderivs.length; i++) {
				densityderivs[i] = new SimpleMatrix(0, 0);
			}

			return densityderivs;
		}

		double[] oldrMags = new double[rarray.length];
		Arrays.fill(oldrMags, 1);

		while (Utils.numNotNull(rarray) > 0) {
			ArrayList<SimpleMatrix> d = new ArrayList<>();

			ArrayList<SimpleMatrix> p = new ArrayList<>();

			for (int i = 0; i < rarray.length; i++) {
				if (rarray[i] != null) {
					d.add(new SimpleMatrix(dirs[i]));
					p.add(D.mult(
							computeResponseVectorsLimited(dirs[i], soln)));
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
				alpha = Utils.filled(solver.numCols(), rhsvec.numCols(), 1);
			}

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					double mag = Utils.mag(rarray[a]);

					for (int i = 0; i < alpha.numRows(); i++) {
						xarray[a] =
								xarray[a].plus(d.get(i).scale(alpha.get(i,
										a)));
						rarray[a] =
								rarray[a].minus(
										p.get(i).scale(alpha.get(i, a)));
					}

					if (mag != mag || oldrMags[a] < mag) {
						throw new IllegalStateException("Thiel has failed!");
					}
					if (mag < 1E5) {
						rarray[a] = null;
					}
					else {
						System.out.println("Thiel convergence test: " + mag);
					}

					oldrMags[a] = mag;
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

			SimpleMatrix beta;
			try {
				beta = solver.solve(rhsvec);
			} catch (SingularMatrixException e) {
				beta = Utils.filled(solver.numCols(), rhsvec.numCols(), 1);
			}

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
		return xarray;
	}


	public static SimpleMatrix xArrayComplementary(SolutionR soln,
												   SimpleMatrix fockderiv) {
		int NOcc = (int) (soln.nElectrons / 2.0);
		if (NOcc == 0) {
			return new SimpleMatrix(0, 0);
		}
		SimpleMatrix x = new SimpleMatrix(NOcc - 1, 1);
		int count1 = 0;

		for (int j = 0; j < NOcc - 1; j++) {
			double element = 0;

			for (int u = 0; u < soln.orbitals.length; u++) {
				for (int v = 0; v < soln.orbitals.length; v++) {
					element += soln.C.get(NOcc - 1, u) * soln.C.get(j, v) *
							fockderiv.get(u, v);
				}
			}
			x.set(count1, 0,
					-element / (soln.E.get(NOcc - 1) - soln.E.get(j)));
			count1++;
		}
		return x;
	}


	public static SimpleMatrix densityDerivativeLimited(SolutionR soln,
														SimpleMatrix x) {

		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix densityMatrixDeriv =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);
		for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
			for (int v = 0; v < densityMatrixDeriv.numCols(); v++) {
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
				densityMatrixDeriv.set(u, v, sum);
			}
		}

		return densityMatrixDeriv;
	}

	public static SimpleMatrix xarrayForIE(SolutionR soln,
										   SimpleMatrix xlimited,
										   SimpleMatrix xcomplementary) {

		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix x = new SimpleMatrix(soln.orbitals.length - 1, 1);

		int count = 0;

		for (int i = xlimited.getNumElements() - NVirt;
			 i < xlimited.getNumElements(); i++) {
			if (i > -1) {
				x.set(count, 0, xlimited.get(i, 0));
				count++;
			}
		}

		for (int i = 0; i < xcomplementary.getNumElements(); i++) {
			x.set(count, 0, xcomplementary.get(i, 0));
			count++;
		}

		return x;
	}


	public static SimpleMatrix HOMOCoefficientDerivativeComplementary(
			SimpleMatrix x, SolutionR soln) {


		SimpleMatrix CDeriv = new SimpleMatrix(1, soln.orbitals.length);

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

			CDeriv.set(0, u, sum);
		}

		return CDeriv;
	}


	public static double MNDOIEDeriv(SolutionR soln, SimpleMatrix coeffDeriv,
									 SimpleMatrix Fderiv) {

		int index = (int) (soln.nElectrons / 2.0) - 1;

		if (index < 0) {
			return 0;
		}

		SimpleMatrix coeff = soln.C.extractVector(true, index);

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
										 SimpleMatrix densityderiv, int Z,
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

		int[][] index = soln.orbsOfAtom;

		SimpleMatrix densityMatrix = soln.densityMatrix();

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
