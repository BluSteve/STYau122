package nddoparam;

import nddoparam.param.ParamDerivative;
import scf.GTO;
import scf.LCGTO;
import scf.OrbitalProperties;
import scf.STO6G;

public class NDDO6G extends STO6G {
	public double zeta, beta, U, p0, p1, p2, D1, D2;
	public double gss, gsp, hsp, gpp, gp2, hp2;
	private NDDOAtom a;
	private NDDO6G[] orbitalArray;

	public NDDO6G(NDDOAtom a, OrbitalProperties orbital, double zeta,
				  double beta, double U) {
		super(zeta, a, orbital);

		this.a = a;
		this.zeta = zeta;
		this.beta = beta;
		this.U = U;
		this.p0 = a.p0;
		this.p1 = a.p1;
		this.p2 = a.p2;
		this.D1 = a.D1;
		this.D2 = a.D2;
		this.gss = a.getParams().getGss();
		this.gsp = a.getParams().getGsp();
		this.hsp = a.getParams().getHsp();
		this.gpp = a.getParams().getGpp();
		this.gp2 = a.getParams().getGp2();
		this.hp2 = 0.5 * (gpp - gp2);
	}

	public NDDO6G(int i, int j, int k, NDDO6G nddo6G) {
		super.i = i;
		super.j = j;
		super.k = k;
		super.L = i + j + k;
		super.coordinates = nddo6G.coordinates;
		this.p0 = nddo6G.p0;
		this.p1 = nddo6G.p1;
		this.p2 = nddo6G.p2;
		this.D1 = nddo6G.D1;
		this.D2 = nddo6G.D2;
	}

	public NDDO6G(NDDO6G nddo6G, double[] coordinates) {
		this(nddo6G.i, nddo6G.j, nddo6G.k, nddo6G);
		super.coordinates = coordinates;
	}

	public static double OneCenterERI(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d) {
		if (a.getL() == b.getL() && a.getL() == 0) {//(ss|??)
			if (c.getL() == d.getL() & c.getL() == 0) {//(gss
				//System.err.println ("gss");
				return a.gss;
			}
			else if ((c.geti() == 1 && d.geti() == 1) ||
					(c.getj() == 1 && d.getj() == 1) ||
					(c.getk() == 1 && d.getk() == 1)) {//gsp
				//System.err.println ("gsp");
				return a.gsp;
			}
			else {
				return 0;
			}
		}
		else if (a.getL() == 1 && b.getL() == 1) {//(pp'|??)

			if ((a.geti() == 1 && b.geti() == 1) ||
					(a.getj() == 1 && b.getj() == 1) ||
					(a.getk() == 1 && b.getk() == 1)) {//(pp|??)
				if (c.getL() == d.getL() && c.getL() == 0) {//gsp
					//System.err.println ("gsp");
					return a.gsp;
				}
				else if ((c.geti() == 1 && d.geti() == 1) ||
						(c.getj() == 1 && d.getj() == 1) ||
						(c.getk() == 1 && d.getk() == 1)) {
					if (a.geti() == c.geti() && a.geti() == 1 ||
							a.getj() == c.getj() && a.getj() == 1 ||
							a.getk() == c.getk() && a.getk() == 1) {//gpp
						//System.err.println ("gpp");
						return a.gpp;
					}
					else {//gpp'
						//System.err.println ("gpp'");
						return a.gp2;
					}
				}

			}
			else if (c.getL() == d.getL() && c.getL() == 1) {//(pp'|p''p''')
				if (Math.abs(LCGTO.getS(a, c) - 1) < 1E-5 &&
						Math.abs(LCGTO.getS(b, d) - 1) < 1E-5) {
					//System.err.println ("hpp");
					return a.hp2;
				}
				else if (Math.abs(LCGTO.getS(a, d) - 1) < 1E-5 &&
						Math.abs(LCGTO.getS(b, c) - 1) < 1E-5) {
					//System.err.println ("hpp");
					return a.hp2;
				}
			}
		}
		else if (a.getL() == 0 && c.getL() == 0 &&
				((b.geti() == 1 && d.geti() == 1) ||
						(b.getj() == 1 && d.getj() == 1) ||
						(b.getk() == 1 && d.getk() == 1))) {
			//System.err.println ("hsp");
			return a.hsp;
		}
		else if (a.getL() == 0 && d.getL() == 0 &&
				((b.geti() == 1 && c.geti() == 1) ||
						(b.getj() == 1 && c.getj() == 1) ||
						(b.getk() == 1 && c.getk() == 1))) {
			//System.err.println ("hsp");
			return a.hsp;
		}
		else if (b.getL() == 0 && c.getL() == 0 &&
				((a.geti() == 1 && d.geti() == 1) ||
						(a.getj() == 1 && d.getj() == 1) ||
						(a.getk() == 1 && d.getk() == 1))) {
			//System.err.println ("hsp");
			return a.hsp;
		}
		else if (b.getL() == 0 && d.getL() == 0 &&
				((a.geti() == 1 && c.geti() == 1) ||
						(a.getj() == 1 && c.getj() == 1) ||
						(a.getk() == 1 && c.getk() == 1))) {
			//System.err.println ("hsp");
			return a.hsp;
		}
		return 0;
	}

	public static double beta(NDDO6G a, NDDO6G b) {
		return 0.5 * (a.beta + b.beta) * LCGTO.getS(a, b);
	}

	public static double betaderiv(NDDO6G a, NDDO6G b, int tau) {
		return 0.5 * (a.beta + b.beta) * LCGTO.getSDeriv(a, b, tau);
	}

	public static double betaparamderiv(NDDO6G a, NDDO6G b, int num,
										int type) {
		return 0.5 * (a.beta + b.beta) * ParamDerivative
				.getSderivfinite(a, b, num, type);
	}

	public static double betaderiv2(NDDO6G a, NDDO6G b, int tau1, int tau2) {
		return 0.5 * (a.beta + b.beta) * LCGTO.getSDeriv2(a, b, tau1, tau2);
	}

	private static double qq(double p01, double p11, double p21, double D11,
							 double D21, double p02, double p12, double p22,
							 double D12, double D22, double R) {
		double a00 = p01 + p02;
		return Math.pow(R * R + a00 * a00, -0.5);
	}

	private static double quz(double p01, double p11, double p21, double D11,
							  double D21, double p02, double p12, double p22,
							  double D12, double D22, double R) {
		double a01 = p01 + p12;
		return 0.5 * Math.pow((R + D12) * (R + D12) + a01 * a01, -0.5)
				- 0.5 * Math.pow((R - D12) * (R - D12) + a01 * a01, -0.5);
	}

	private static double qQpipi(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R) {
		double a02 = p01 + p22;
		return 0.5 * Math.pow(R * R + 4 * D22 * D22 + a02 * a02, -0.5)
				- 0.5 * Math.pow(R * R + a02 * a02, -0.5);
	}

	private static double qQzz(double p01, double p11, double p21, double D11,
							   double D21, double p02, double p12, double p22,
							   double D12, double D22, double R) {
		double a02 = p01 + p22;
		return 0.25 * Math.pow((R + 2 * D22) * (R + 2 * D22) + a02 * a02, -0.5)
				+
				0.25 * Math.pow((R - 2 * D22) * (R - 2 * D22) + a02 * a02,
						-0.5)
				- 0.5 * Math.pow(R * R + a02 * a02, -0.5);
	}

	private static double upiupi(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R) {
		double a11 = p11 + p12;
		return 0.5 *
				Math.pow(R * R + (D11 - D12) * (D11 - D12) + a11 * a11, -0.5)
				- 0.5 *
				Math.pow(R * R + (D11 + D12) * (D11 + D12) + a11 * a11, -0.5);
	}

	private static double uzuz(double p01, double p11, double p21, double D11,
							   double D21, double p02, double p12, double p22,
							   double D12, double D22, double R) {
		double a11 = p11 + p12;
		return 0.25 *
				Math.pow((R + D11 - D12) * (R + D11 - D12) + a11 * a11, -0.5)
				- 0.25 *
				Math.pow((R + D11 + D12) * (R + D11 + D12) + a11 * a11, -0.5)
				- 0.25 *
				Math.pow((R - D11 - D12) * (R - D11 - D12) + a11 * a11, -0.5)
				+ 0.25 *
				Math.pow((R - D11 + D12) * (R - D11 + D12) + a11 * a11, -0.5);
	}

	private static double upiQpiz(double p01, double p11, double p21,
								  double D11, double D21, double p02,
								  double p12, double p22, double D12,
								  double D22, double R) {
		double a12 = p11 + p22;
		return -0.25 * Math.pow(
				(R - D22) * (R - D22) + (D11 - D22) * (D11 - D22) + a12 * a12,
				-0.5)
				+ 0.25 * Math.pow(
				(R - D22) * (R - D22) + (D11 + D22) * (D11 + D22) + a12 * a12,
				-0.5)
				+ 0.25 * Math.pow(
				(R + D22) * (R + D22) + (D11 - D22) * (D11 - D22) + a12 * a12,
				-0.5)
				- 0.25 * Math.pow(
				(R + D22) * (R + D22) + (D11 + D22) * (D11 + D22) + a12 * a12,
				-0.5);
	}

	private static double uzQpipi(double p01, double p11, double p21,
								  double D11, double D21, double p02,
								  double p12, double p22, double D12,
								  double D22, double R) {
		double a12 = p11 + p22;
		return -0.25 *
				Math.pow((R + D11) * (R + D11) + 4 * D22 * D22 + a12 * a12,
						-0.5)
				+ 0.25 *
				Math.pow((R - D11) * (R - D11) + 4 * D22 * D22 + a12 * a12,
						-0.5)
				+ 0.25 * Math.pow((R + D11) * (R + D11) + a12 * a12, -0.5)
				- 0.25 * Math.pow((R - D11) * (R - D11) + a12 * a12, -0.5);
	}

	private static double uzQzz(double p01, double p11, double p21, double D11,
								double D21, double p02, double p12, double p22,
								double D12, double D22, double R) {
		double a12 = p11 + p22;
		return -0.125 *
				Math.pow((R + D11 - 2 * D22) * (R + D11 - 2 * D22) + a12 * a12,
						-0.5)
				+ 0.125 *
				Math.pow((R - D11 - 2 * D22) * (R - D11 - 2 * D22) + a12 * a12,
						-0.5)
				- 0.125 *
				Math.pow((R + D11 + 2 * D22) * (R + D11 + 2 * D22) + a12 * a12,
						-0.5)
				+ 0.125 *
				Math.pow((R - D11 + 2 * D22) * (R - D11 + 2 * D22) + a12 * a12,
						-0.5)
				+ 0.25 * Math.pow((R + D11) * (R + D11) + a12 * a12, -0.5)
				- 0.25 * Math.pow((R - D11) * (R - D11) + a12 * a12, -0.5);
	}

	private static double QpipiQpipi(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double R) {
		double a22 = p21 + p22;
		return 0.125 *
				Math.pow(R * R + 4 * (D21 - D22) * (D21 - D22) + a22 * a22,
						-0.5)
				+ 0.125 *
				Math.pow(R * R + 4 * (D21 + D22) * (D21 + D22) + a22 * a22,
						-0.5)
				- 0.25 * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -0.5)
				- 0.25 * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -0.5)
				+ 0.25 * Math.pow(R * R + a22 * a22, -0.5);
	}

	private static double QxxQyy(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R) {
		double a22 = p21 + p22;
		return 0.25 *
				Math.pow(R * R + 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22,
						-0.5)
				- 0.25 * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -0.5)
				- 0.25 * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -0.5)
				+ 0.25 * Math.pow(R * R + a22 * a22, -0.5);
	}

	private static double QpipiQzz(double p01, double p11, double p21,
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R) {
		double a22 = p21 + p22;
		return 0.125 * Math.pow(
				(R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 + a22 * a22,
				-0.5)
				+ 0.125 * Math.pow(
				(R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 + a22 * a22,
				-0.5)
				- 0.125 *
				Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -0.5)
				- 0.125 *
				Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -0.5)
				- 0.25 * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -0.5)
				+ 0.25 * Math.pow(R * R + a22 * a22, -0.5);
	}

	private static double QzzQzz(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R) {
		double a22 = p21 + p22;
		return 0.0625 * Math.pow(
				(R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) + a22 * a22,
				-0.5)
				+ 0.0625 * Math.pow(
				(R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) + a22 * a22,
				-0.5)
				+ 0.0625 * Math.pow(
				(R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) + a22 * a22,
				-0.5)
				+ 0.0625 * Math.pow(
				(R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) + a22 * a22,
				-0.5)
				- 0.125 *
				Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -0.5)
				- 0.125 *
				Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -0.5)
				- 0.125 *
				Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -0.5)
				- 0.125 *
				Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -0.5)
				+ 0.25 * Math.pow(R * R + a22 * a22, -0.5);
	}

	private static double QpizQpiz(double p01, double p11, double p21,
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R) {
		double a22 = p21 + p22;
		return 0.125 * Math.pow(
				(R + D21 - D22) * (R + D21 - D22) + (D21 - D22) * (D21 - D22) +
						a22 * a22, -0.5)
				- 0.125 * Math.pow(
				(R + D21 - D22) * (R + D21 - D22) + (D21 + D22) * (D21 + D22) +
						a22 * a22, -0.5)
				- 0.125 * Math.pow(
				(R + D21 + D22) * (R + D21 + D22) + (D21 - D22) * (D21 - D22) +
						a22 * a22, -0.5)
				+ 0.125 * Math.pow(
				(R + D21 + D22) * (R + D21 + D22) + (D21 + D22) * (D21 + D22) +
						a22 * a22, -0.5)
				- 0.125 * Math.pow(
				(R - D21 - D22) * (R - D21 - D22) + (D21 - D22) * (D21 - D22) +
						a22 * a22, -0.5)
				+ 0.125 * Math.pow(
				(R - D21 - D22) * (R - D21 - D22) + (D21 + D22) * (D21 + D22) +
						a22 * a22, -0.5)
				+ 0.125 * Math.pow(
				(R - D21 + D22) * (R - D21 + D22) + (D21 - D22) * (D21 - D22) +
						a22 * a22, -0.5)
				- 0.125 * Math.pow(
				(R - D21 + D22) * (R - D21 + D22) + (D21 + D22) * (D21 + D22) +
						a22 * a22, -0.5);
	}

	private static double ssss(double p01, double p11, double p21, double D11,
							   double D21, double p02, double p12, double p22,
							   double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double ssppippi(double p01, double p11, double p21,
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double sspzpz(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double ppippiss(double p01, double p11, double p21,
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	private static double pzpzss(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	private static double ppippippippi(double p01, double p11, double p21,
									   double D11, double D21, double p02,
									   double p12, double p22, double D12,
									   double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QpipiQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						R);
	}

	private static double pxpxpypy(double p01, double p11, double p21,
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QxxQyy(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double ppippipzpz(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QpipiQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double pzpzppippi(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QpipiQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	private static double pzpzpzpz(double p01, double p11, double p21,
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R) {
		return qq(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) +
				qQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				QzzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double spzss(double p01, double p11, double p21, double D11,
								double D21, double p02, double p12, double p22,
								double D12, double D22, double R) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	private static double spzppippi(double p01, double p11, double p21,
									double D11, double D21, double p02,
									double p12, double p22, double D12,
									double D22, double R) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				uzQpipi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double spzpzpz(double p01, double p11, double p21,
								  double D11, double D21, double p02,
								  double p12, double p22, double D12,
								  double D22, double R) {
		return -quz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R) +
				uzQzz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double ssspz(double p01, double p11, double p21, double D11,
								double D21, double p02, double p12, double p22,
								double D12, double D22, double R) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double ppippispz(double p01, double p11, double p21,
									double D11, double D21, double p02,
									double p12, double p22, double D12,
									double D22, double R) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) -
				uzQpipi(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	private static double pzpzspz(double p01, double p11, double p21,
								  double D11, double D21, double p02,
								  double p12, double p22, double D12,
								  double D22, double R) {
		return quz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R) -
				uzQzz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	private static double sppisppi(double p01, double p11, double p21,
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R) {
		return upiupi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double spzspz(double p01, double p11, double p21,
								 double D11,
								 double D21, double p02, double p12,
								 double p22,
								 double D12, double D22, double R) {
		return uzuz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double sppippipz(double p01, double p11, double p21,
									double D11, double D21, double p02,
									double p12, double p22, double D12,
									double D22, double R) {
		return upiQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double ppipzsppi(double p01, double p11, double p21,
									double D11, double D21, double p02,
									double p12, double p22, double D12,
									double D22, double R) {
		return -upiQpiz(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R);
	}

	private static double ppipzppipz(double p01, double p11, double p21,
									 double D11, double D21, double p02,
									 double p12, double p22, double D12,
									 double D22, double R) {
		return QpizQpiz(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R);
	}

	private static double pxpypxpy(double p01, double p11, double p21,
								   double D11, double D21, double p02,
								   double p12, double p22, double D12,
								   double D22, double R) {
		return 0.5 *
				(ppippippippi(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22,
						R) -
						pxpxpypy(p01, p11, p21, D11, D21, p02, p12, p22, D12,
								D22, R));
	}

	public static double LocalTwoCenterERI(NDDO6G a, NDDO6G b, NDDO6G c,
										   NDDO6G d) {

		double R = GTO.R(a.getCoords(), c.getCoords());
		//(??|??)
		switch (a.L) {

			case 0://(s?|??)

				switch (b.L) {

					case 0: //(ss|??)

						switch (c.L) {

							case 0: //(ss|s?);

								switch (d.L) {

									case 0://(ss|ss)
										return ssss(a.p0, a.p1, a.p2, a.D1,
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2, R);

									case 1:
										if (d.k == 1) {//(ss|spz)
											return ssspz(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2, R);
										}
										else {//(ss|sppi) = 0
											return 0;
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							case 1: //(ss|p?)
								if (c.k == 1) {//(ss|pz?)

									switch (d.L) {

										case 0://(ss|pzs)
											return ssspz(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2, R);

										case 1:
											if (d.k == 1) {//(ss|pzpz)
												return sspzpz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);
											}
											else {//(ss|pzppi) = 0
												return 0;
											}
										default:
											return 0;
									}
								}
								else {//(ss|ppi?)

									if (d.L == 1 && d.k == 0 && c.i == d.i &&
											c.j == d.j) {//(ss|ppippi)
										return ssppippi(a.p0, a.p1, a.p2, a.D1,
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2, R);
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

						if (b.k == 1) {//(spz|??)

							switch (c.L) {

								case 0://(spz|s?)

									switch (d.L) {

										case 0://(spz|ss)
											return spzss(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2, R);

										case 1:
											if (d.k == 1) {//(spz|spz)
												return spzspz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.k == 1) {//(spz|pz?)

										switch (d.L) {

											case 0://(spz|pzs)
												return spzspz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);

											case 1:
												if (d.k == 1) {//(spz|pzpz)
													return spzpzpz(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {//(spz|pzppi) = 0
													return 0;
												}
										}
									}
									else {//(spz|ppi?)
										if (d.i == c.i && d.j == c.j &&
												d.k == 0) {
											return spzppippi(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
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

							switch (c.L) {
								case 0://(sppi|s?)
									if (d.i == b.i && d.j == b.j &&
											d.k == 0) {//(sppi|sppi)
										return sppisppi(a.p0, a.p1, a.p2, a.D1,
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2, R);
									}
									else {
										return 0;
									}
								case 1:
									if (c.k == 1) {
										if (d.i == b.i && d.j == b.j &&
												d.k == 0) {//(sppi|pzppi)
											return sppippipz(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									}
									else {
										if (c.i == b.i && c.j == b.j &&
												c.k == 0) {//(sppi|ppi?)
											switch (d.L) {
												case 0:
													return sppisppi(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												case 1:
													if (d.k == 1) {
														return sppippipz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
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
				if (a.k == 1) {//(pz?|??)
					switch (b.L) {
						case 0:
							switch (c.L) {

								case 0://(pzs|s?)

									switch (d.L) {

										case 0://(pzs|ss)
											return spzss(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2, R);

										case 1:
											if (d.k == 1) {//(pzs|spz)
												return spzspz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.k == 1) {//(pzs|pz?)

										switch (d.L) {

											case 0://(pzs|pzs)
												return spzspz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);

											case 1:
												if (d.k == 1) {//(pzs|pzpz)
													return spzpzpz(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {//(pzs|pzppi) = 0
													return 0;
												}
										}
									}
									else {//(pzs|ppi?)
										if (d.i == c.i && d.j == c.j &&
												d.k == 0) {
											return spzppippi(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
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

							if (b.k == 1) {//(pzpz|??)

								switch (c.L) {

									case 0://(pzpz|s?)

										switch (d.L) {

											case 0://(pzpz|ss)
												return pzpzss(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);

											case 1:
												if (d.k == 1) {//(pzpz|spz)
													return pzpzspz(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {
													return 0;
												}
										}

									case 1:
										if (c.k == 1) {//(pzpz|pz?)

											switch (d.L) {

												case 0://(pzpz|pzs)
													return pzpzspz(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);

												case 1:
													if (d.k == 1) {//(pzpz
														// |pzpz)
														return pzpzpzpz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
													else {//(pzpz|pzppi) = 0
														return 0;
													}
											}
										}
										else {//(pzpz|ppi?)
											if (d.i == c.i && d.j == c.j &&
													d.k == 0) {
												return pzpzppippi(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R);
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

								switch (c.L) {
									case 0://(pzppi|s?)
										if (d.i == b.i && d.j == b.j &&
												d.k == 0) {//(pzppi|sppi)
											return ppipzsppi(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									case 1:
										if (c.k == 1) {
											if (d.i == b.i && d.j == b.j &&
													d.k == 0) {//(pzppi|pzppi)
												return ppipzppipz(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R);
											}
											else {
												return 0;
											}
										}
										else {
											if (c.i == b.i && c.j == b.j &&
													c.k == 0) {//(pzppi|ppi?)
												switch (d.L) {
													case 0:
														return ppipzsppi(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													case 1:
														if (d.k == 1) {
															return ppipzppipz(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	R);
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

					switch (b.L) {
						case 0://(ppis|??)

							switch (c.L) {
								case 0://(ppis|s?)
									if (d.i == a.i && d.j == a.j &&
											d.k == 0) {//(ppis|sppi)
										return sppisppi(a.p0, a.p1, a.p2, a.D1,
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2, R);
									}
									else {
										return 0;
									}
								case 1:
									if (c.k == 1) {
										if (d.i == a.i && d.j == a.j &&
												d.k == 0) {//(ppis|pzppi)
											return sppippipz(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									}
									else {
										if (c.i == a.i && c.j == a.j &&
												c.k == 0) {//(ppis|ppi?)
											switch (d.L) {
												case 0:
													return sppisppi(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												case 1:
													if (d.k == 1) {
														return sppippipz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
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
							if (b.k == 1) {//(ppipz|??)
								switch (c.L) {
									case 0://(ppipz|s?)
										if (d.i == a.i && d.j == a.j &&
												d.k == 0) {//(ppipz|sppi)
											return ppipzsppi(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									case 1:
										if (c.k == 1) {
											if (d.i == a.i && d.j == a.j &&
													d.k == 0) {//(ppipz|pzppi)
												return ppipzppipz(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R);
											}
											else {
												return 0;
											}
										}
										else {
											if (c.i == a.i && c.j == a.j &&
													c.k == 0) {//(ppipz|ppi?)
												switch (d.L) {
													case 0:
														return ppipzsppi(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													case 1:
														if (d.k == 1) {
															return ppipzppipz(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	R);
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

								switch (c.L) {
									case 0://(ppippi|s?)
										switch (d.L) {
											case 0://(ppippi|ss)
												if (a.i == b.i && a.j == b.j) {
													return ppippiss(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {
													return 0;
												}
											case 1:
												if (d.k == 1 && a.i == b.i &&
														a.j == b.j) {
													return ppippispz(a.p0,
															a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {
													return 0;
												}
										}

									case 1:
										if (c.k == 1) {
											switch (d.L) {
												case 0://(ppippi|pzs)
													if (a.i == b.i &&
															a.j == b.j) {
														return ppippispz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
													else {
														return 0;
													}

												case 1:
													if (d.k == 1 &&
															a.i == b.i && a.j ==
															b.j) {//(ppippi
														// |pzpz)
														return ppippipzpz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
													else {
														return 0;
													}
											}
										}
										else {
											if (a.i == b.i && a.j ==
													b.j) {//(pxpx|??) or
												// (pypy|??)

												if (c.L == d.L && c.i == d.i &&
														c.j == d.j &&
														c.k == 0) {
													if (a.i == c.i) {
														return ppippippippi(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, R);
													}
													else {
														return pxpxpypy(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
												}
												else {
													return 0;
												}

											}
											else {//(pxpy|??) or (pypx|??)
												if (c.L == d.L && c.i != d.i &&
														c.j != d.j &&
														c.k == 0) {
													return pxpypxpy(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
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

	public static double getG(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d) {
		double[] coeffA = a.decomposition(a.coordinates, c.coordinates);
		double[] coeffB = b.decomposition(a.coordinates, c.coordinates);
		double[] coeffC = c.decomposition(a.coordinates, c.coordinates);
		double[] coeffD = d.decomposition(a.coordinates, c.coordinates);

		double sum2 = 0;

		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {

						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] !=
								0) {
							sum2 += coeffA[i] * coeffB[j] * coeffC[k] *
									coeffD[l] *
									LocalTwoCenterERI(A[i], B[j], C[k], D[l]) *
									27.21;
						}


					}
				}
			}
		}


		return sum2;
	}

	public NDDOAtom getAtom() {
		return this.a;
	}

	public double U() {
		return this.U;
	}

	public double beta() {
		return this.beta;
	}

	public double[] decomposition(double[] point1, double[] point2) {

		if (this.L == 0) {
			return new double[]{1};
		}
		else if (this.L == 1) {

			double[] zloc = new double[3];
			double val = GTO.R(point1, point2);
			for (int i = 0; i < 3; i++) {
				zloc[i] = (point2[i] - point1[i]) / val;
			}

			double[] yloc = new double[3];
			double scale = Math.sqrt(zloc[0] * zloc[0] + zloc[1] * zloc[1]);

			yloc[0] = -zloc[1] / scale;
			yloc[1] = zloc[0] / scale;
			yloc[2] = 0;

			if (scale == 0) {
				yloc[0] = 0;
				yloc[1] = 1;
			}

			double[] xloc = new double[]{-yloc[1] * zloc[2], yloc[0] * zloc[2],
					zloc[0] * yloc[1] - zloc[1] * yloc[0]};

			if (this.i == 1) {
				return new double[]{xloc[0], yloc[0], zloc[0]};
			}
			else if (this.j == 1) {
				return new double[]{xloc[1], yloc[1], zloc[1]};
			}
			else if (this.k == 1) {
				return new double[]{xloc[2], yloc[2], zloc[2]};
			}
		}
		return null;
	}

	public double[] decomposition2(double[] point1, double[] point2) {

		if (this.L == 0) {
			return new double[]{1};
		}

		double x = point2[0] - point1[0];

		double y = point2[1] - point1[1];

		double z = point2[2] - point1[2];

		double Rxz = Math.sqrt(x * x + z * z);

		double R = GTO.R(point1, point2);

		if (this.L == 1) {
			if (this.i == 1) {
				if (Rxz == 0) {
					return new double[]{1, 0, 0};
				}
				return new double[]{x * y / (R * Rxz), -z / Rxz, x / R};
			}
			else if (this.j == 1) {
				return new double[]{-Rxz / R, 0, y / R};
			}
			else {
				if (Rxz == 0) {
					return new double[]{0, 1, 0};
				}
				return new double[]{y * z / (R * Rxz), x / Rxz, z / R};
			}
		}

		return null;


	}

	public NDDO6G[] orbitalArray() {
		if (this.orbitalArray == null) {
			if (this.L == 0) {
				orbitalArray = new NDDO6G[]{new NDDO6G(0, 0, 0, this)};
			}
			else if (this.L == 1) {
				orbitalArray = new NDDO6G[]{
						new NDDO6G(1, 0, 0, this),
						new NDDO6G(0, 1, 0, this),
						new NDDO6G(0, 0, 1, this)};
			}
		}

		return this.orbitalArray;
	}
}
