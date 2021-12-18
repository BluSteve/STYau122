package nddo.math;

public class Multipoles {

	public static double qq(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return Math.pow(R * R + (p01 + p02) * (p01 + p02), -0.5);
	}

	public static double f0 (int a2, int b2, double Dn2, double p01, double pn2, double R) {
		return Math.pow((R + a2 * Dn2) * (R + a2 * Dn2) + b2 * b2 * Dn2 * Dn2 + (p01 + pn2) * (p01 + pn2), -0.5);
	}

	public static double quz(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.5 * (
				+ f0(+1, 0, D12, p01, p12, R)
						- f0(-1, 0, D12, p01, p12, R)
		);
	}

	public static double qQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.5 * (
				+ f0(0, 2, D22, p01, p22, R)
						- f0(0, 0, D22, p01, p22, R)
		);
	}

	public static double qQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.25 * (
				+ f0(+2, 0, D22, p01, p22, R)
						+ f0(-2, 0, D22, p01, p22, R)
						- 2 * f0(0, 0, D22, p01, p22, R)
		);
	}

	public static double f1 (int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R) {
		return Math.pow((R + a1 * Dn1 + a2 * Dm2) * (R + a1 * Dn1 + a2 * Dm2) + (b1 * Dn1 + b2 * Dm2) * (b1 * Dn1 + b2 * Dm2) + (pn1 + pm2) * (pn1 + pm2), -0.5);
	}

	public static double f2 (int a1, int a2, int b1, int b2, double D21, double D22, double p21, double p22, double R) {
		return Math.pow((R + a1 * D21 + a2 * D22) * (R + a1 * D21 + a2 * D22) + (b1 * D21) * (b1 * D21) + (b2 * D22) * (b2 * D22) + (p21 + p22) * (p21 + p22), -0.5);
	}

	public static double upiupi(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.5 * (
				+ f1(0, 0, 1, -1, D11, D12, p11, p12, R)
						- f1(0, 0, 1, +1, D11, D12, p11, p12, R)
		);
	}

	public static double uzuz(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.25 * (
				+ f1(+1, -1, 0, 0, D11, D12, p11, p12, R)
						+ f1(-1, +1, 0, 0, D11, D12, p11, p12, R)
						- f1(-1, -1, 0, 0, D11, D12, p11, p12, R)
						- f1(+1, +1, 0, 0, D11, D12, p11, p12, R)
		);
	}

	public static double upiQpiz(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.25 * (
				+ f1(0, -1, +1, +1, D11, D22, p11, p22, R)
						- f1(0, -1, +1, -1, D11, D22, p11, p22, R)
						+ f1(0, +1, +1, -1, D11, D22, p11, p22, R)
						- f1(0, +1, +1, +1, D11, D22, p11, p22, R)
		);
	}

	public static double uzQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.25 * (
				- f1(+1, 0, 0, 2, D11, D22, p11, p22, R)
						+ f1(-1, 0, 0, 2, D11, D22, p11, p22, R)
						+ f1(+1, 0, 0, 0, D11, D22, p11, p22, R)
						- f1(-1, 0, 0, 0, D11, D22, p11, p22, R)
		);
	}

	public static double uzQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.125 * (
				- f1(+1, -2, 0, 0, D11, D22, p11, p22, R)
						+ f1(-1, -2, 0, 0, D11, D22, p11, p22, R)
						- f1(+1, +2, 0, 0, D11, D22, p11, p22, R)
						+ f1(-1, +2, 0, 0, D11, D22, p11, p22, R)
						+ 2 * f1(+1, 0, 0, 0, D11, D22, p11, p22, R)
						- 2 * f1(-1, 0, 0, 0, D11, D22, p11, p22, R)
		);
	}

	public static double QpipiQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.125 * (
				+ f1(0, 0, +2, -2, D21, D22, p21, p22, R)
						+ f1(0, 0, +2, +2, D21, D22, p21, p22, R)
						- 2 * f1(0, 0, +2, 0, D21, D22, p21, p22, R)
						- 2 * f1(0, 0, 0, +2, D21, D22, p21, p22, R)
						+ 2 * f1(0, 0, 0, 0, D21, D22, p21, p22, R)
		);
	}

	public static double QxxQyy(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.25 * (
				+ f2(0, 0, 2, 2, D21, D22, p21, p22, R)
						- f2(0, 0, 2, 0, D21, D22, p21, p22, R)
						- f2(0, 0, 0, 2, D21, D22, p21, p22, R)
						+ f2(0, 0, 0, 0, D21, D22, p21, p22, R)
		);
	}

	public static double QpipiQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.125 * (
				+ f1(0, -2, +2, 0, D21, D22, p21, p22, R)
						+ f1(0, +2, +2, 0, D21, D22, p21, p22, R)
						- f1(0, -2, 0, 0, D21, D22, p21, p22, R)
						- f1(0, +2, 0, 0, D21, D22, p21, p22, R)
						- 2 * f1(0, 0, 2, 0, D21, D22, p21, p22, R)
						+ 2 * f1(0, 0, 0, 0, D21, D22, p21, p22, R)
		);
	}

	public static double QzzQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.0625 * (
				+ f1(+2, -2, 0, 0, D21, D22, p21, p22, R)
						+ f1(+2, +2, 0, 0, D21, D22, p21, p22, R)
						+ f1(-2, -2, 0, 0, D21, D22, p21, p22, R)
						+ f1(-2, +2, 0, 0, D21, D22, p21, p22, R)
						- 2 * f1(+2, 0, 0, 0, D21, D22, p21, p22, R)
						- 2 * f1(-2, 0, 0, 0, D21, D22, p21, p22, R)
						- 2 * f1(0, +2, 0, 0, D21, D22, p21, p22, R)
						- 2 * f1(0, -2, 0, 0, D21, D22, p21, p22, R)
						+ 4 * f1(0, 0, 0, 0, D21, D22, p21, p22, R)
		);
	}

	public static double QpizQpiz(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R) {
		return 0.125 * (
				+ f1(+1, -1, +1, -1, D21, D22, p21, p22, R)
						- f1(+1, -1, +1, +1, D21, D22, p21, p22, R)
						- f1(+1, +1, +1, -1, D21, D22, p21, p22, R)
						+ f1(+1, +1, +1, +1, D21, D22, p21, p22, R)
						- f1(-1, -1, +1, -1, D21, D22, p21, p22, R)
						+ f1(-1, -1, +1, +1, D21, D22, p21, p22, R)
						+ f1(-1, +1, +1, -1, D21, D22, p21, p22, R)
						- f1(-1, +1, +1, +1, D21, D22, p21, p22, R)
		);
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
	public static double qqgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {

		return (B[tau] - A[tau]) * Math.pow(R * R + (p01 + p02) * (p01 + p02), -1.5);
	}

	public static double f0gd (int a2, int b2, double Dn2, double p01, double pn2, double R, double[] A, double[] B, int tau) {
		return (B[tau] - A[tau]) / R * (R + a2 * Dn2) * Math.pow(f0(a2, b2, Dn2, p01, pn2, R), 3);
	}

	public static double quzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.5 * (
				+ f0gd(+1, 0, D12, p01, p12, R, A, B, tau)
						- f0gd(-1, 0, D12, p01, p12, R, A, B, tau)
		);
	}

	public static double qQpipigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.5 * (
				+ f0gd(0, 2, D22, p01, p22, R, A, B, tau)
						- f0gd(0, 0, D22, p01, p22, R, A, B, tau)
		);
	}

	public static double qQzzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (
				+ f0gd(+2, 0, D22, p01, p22, R, A, B, tau)
						+ f0gd(-2, 0, D22, p01, p22, R, A, B, tau)
						- 2 * f0gd(0, 0, D22, p01, p22, R, A, B, tau)
		);
	}

	public static double f1gd (int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R, double[] A, double[] B, int tau) {
		return (B[tau] - A[tau]) / R * (R + a1 * Dn1 + a2 * Dm2) * Math.pow(f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
	}

	public static double f2gd (int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R, double[] A, double[] B, int tau) {
		return (B[tau] - A[tau]) / R * (R + a1 * Dn1 + a2 * Dm2) * Math.pow(f2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
	}

	public static double upiupigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.5 * (
				+ f1gd(0, 0, 1, -1, D11, D12, p11, p12, R, A, B, tau)
						- f1gd(0, 0, 1, +1, D11, D12, p11, p12, R, A, B, tau)
		);
	}

	public static double uzuzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (
				+ f1gd(+1, -1, 0, 0, D11, D12, p11, p12, R, A, B, tau)
						+ f1gd(-1, +1, 0, 0, D11, D12, p11, p12, R, A, B, tau)
						- f1gd(-1, -1, 0, 0, D11, D12, p11, p12, R, A, B, tau)
						- f1gd(+1, +1, 0, 0, D11, D12, p11, p12, R, A, B, tau)
		);
	}

	public static double upiQpizgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (
				+ f1gd(0, -1, +1, +1, D11, D22, p11, p22, R, A, B, tau)
						- f1gd(0, -1, +1, -1, D11, D22, p11, p22, R, A, B, tau)
						+ f1gd(0, +1, +1, -1, D11, D22, p11, p22, R, A, B, tau)
						- f1gd(0, +1, +1, +1, D11, D22, p11, p22, R, A, B, tau)
		);
	}

	public static double uzQpipigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (
				- f1gd(+1, 0, 0, 2, D11, D22, p11, p22, R, A, B, tau)
						+ f1gd(-1, 0, 0, 2, D11, D22, p11, p22, R, A, B, tau)
						+ f1gd(+1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau)
						- f1gd(-1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau)
		);
	}

	public static double uzQzzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.125 * (
				- f1gd(+1, -2, 0, 0, D11, D22, p11, p22, R, A, B, tau)
						+ f1gd(-1, -2, 0, 0, D11, D22, p11, p22, R, A, B, tau)
						- f1gd(+1, +2, 0, 0, D11, D22, p11, p22, R, A, B, tau)
						+ f1gd(-1, +2, 0, 0, D11, D22, p11, p22, R, A, B, tau)
						+ 2 * f1gd(+1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau)
						- 2 * f1gd(-1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau)
		);
	}

	public static double QpipiQpipigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.125 * (
				+ f1gd(0, 0, +2, -2, D21, D22, p21, p22, R, A, B, tau)
						+ f1gd(0, 0, +2, +2, D21, D22, p21, p22, R, A, B, tau)
						- 2 * f1gd(0, 0, +2, 0, D21, D22, p21, p22, R, A, B, tau)
						- 2 * f1gd(0, 0, 0, +2, D21, D22, p21, p22, R, A, B, tau)
						+ 2 * f1gd(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau)
		);
	}

	public static double QxxQyygd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (
				+ f2gd(0, 0, 2, 2, D21, D22, p21, p22, R, A, B, tau)
						- f2gd(0, 0, 2, 0, D21, D22, p21, p22, R, A, B, tau)
						- f2gd(0, 0, 0, 2, D21, D22, p21, p22, R, A, B, tau)
						+ f2gd(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau)
		);
	}

	public static double QpipiQzzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.125 * (
				+ f1gd(0, -2, +2, 0, D21, D22, p21, p22, R, A, B, tau)
						+ f1gd(0, +2, +2, 0, D21, D22, p21, p22, R, A, B, tau)
						- f1gd(0, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						- f1gd(0, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						- 2 * f1gd(0, 0, 2, 0, D21, D22, p21, p22, R, A, B, tau)
						+ 2 * f1gd(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau)
		);
	}

	public static double QzzQzzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.0625 * (
				+ f1gd(+2, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						+ f1gd(+2, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						+ f1gd(-2, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						+ f1gd(-2, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						- 2 * f1gd(+2, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						- 2 * f1gd(-2, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						- 2 * f1gd(0, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						- 2 * f1gd(0, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau)
						+ 4 * f1gd(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau)
		);
	}

	public static double QpizQpizgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.125 * (
				+ f1gd(+1, -1, +1, -1, D21, D22, p21, p22, R, A, B, tau)
						- f1gd(+1, -1, +1, +1, D21, D22, p21, p22, R, A, B, tau)
						- f1gd(+1, +1, +1, -1, D21, D22, p21, p22, R, A, B, tau)
						+ f1gd(+1, +1, +1, +1, D21, D22, p21, p22, R, A, B, tau)
						- f1gd(-1, -1, +1, -1, D21, D22, p21, p22, R, A, B, tau)
						+ f1gd(-1, -1, +1, +1, D21, D22, p21, p22, R, A, B, tau)
						+ f1gd(-1, +1, +1, -1, D21, D22, p21, p22, R, A, B, tau)
						- f1gd(-1, +1, +1, +1, D21, D22, p21, p22, R, A, B, tau)
		);
	}

	public static double ssssgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ssppippigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double sspzpzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ppippissgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double pzpzssgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double ppippippippigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) + QpipiQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double pxpxpypygd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) + QxxQyygd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ppippipzpzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) + QpipiQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double pzpzppippigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) + QpipiQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double pzpzpzpzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) + qQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) + QzzQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double spzssgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return -quzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double spzppippigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return -quzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) + uzQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double spzpzpzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return -quzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) + uzQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ssspzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return quzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ppippispzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return quzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) - uzQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double pzpzspzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return quzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) - uzQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double sppisppigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return upiupigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double spzspzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return uzuzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double sppippipzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return upiQpizgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ppipzsppigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return -upiQpizgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double ppipzppipzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return QpizQpizgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double pxpypxpygd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.5 * (ppippippippigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) - pxpxpypygd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau));
	}

	public static double qqg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		double sum = 0;

		if (tau1 == tau2) {

			sum = - Math.pow(R * R + (p01 + p02) * (p01 + p02), -1.5);
		}
		return sum + 3 * (B[tau1] - A[tau1]) * (B[tau2] - A[tau2]) * Math.pow(R * R + (p01 + p02) * (p01 + p02), -2.5);
	}

	public static double f0g2d (int a2, int b2, double Dn2, double p01, double pn2, double R, double[] A, double[] B, int tau1, int tau2) {

		double f0 = f0(a2, b2, Dn2, p01, pn2, R);

		double num = (B[tau1] - A[tau1]) * (B[tau2] - A[tau2]) / (R * R) * Math.pow(f0, 3) * (3 * Math.pow(R + a2 * Dn2, 2) * f0 * f0 + a2 * Dn2 / R);

		if (tau1 == tau2) {
			num += - (R + a2 * Dn2) * Math.pow(f0, 3) / R;
		}

		return num;
	}

	public static double quzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.5 * (
				+ f0g2d(+1, 0, D12, p01, p12, R, A, B, tau1, tau2)
						- f0g2d(-1, 0, D12, p01, p12, R, A, B, tau1, tau2)
		);
	}

	public static double qQpipig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.5 * (
				+ f0g2d(0, 2, D22, p01, p22, R, A, B, tau1, tau2)
						- f0g2d(0, 0, D22, p01, p22, R, A, B, tau1, tau2)
		);
	}

	public static double qQzzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.25 * (
				+ f0g2d(+2, 0, D22, p01, p22, R, A, B, tau1, tau2)
						+ f0g2d(-2, 0, D22, p01, p22, R, A, B, tau1, tau2)
						- 2 * f0g2d(0, 0, D22, p01, p22, R, A, B, tau1, tau2)
		);
	}

	public static double f1g2d (int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R, double[] A, double[] B, int tau1, int tau2) {

		double f1 = f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);

		double num = (B[tau1] - A[tau1]) * (B[tau2] - A[tau2]) / (R * R) * Math.pow(f1, 3) * (3 * Math.pow(R + a1 * Dn1 + a2 * Dm2, 2) * f1 * f1 + (a1 * Dn1 + a2 * Dm2) / R);

		if (tau1 == tau2) {
			num += - (R + a1 * Dn1 + a2 * Dm2) * Math.pow(f1, 3) / R;
		}

		return num;
	}

	public static double f2g2d (int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R, double[] A, double[] B, int tau1, int tau2) {

		double f2 = f2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);

		double num = (B[tau1] - A[tau1]) * (B[tau2] - A[tau2]) / (R * R) * Math.pow(f2, 3) * (3 * Math.pow(R + a1 * Dn1 + a2 * Dm2, 2) * f2 * f2 + (a1 * Dn1 + a2 * Dm2) / R);

		if (tau1 == tau2) {
			num += - (R + a1 * Dn1 + a2 * Dm2) * Math.pow(f2, 3) / R;
		}

		return num;
	}

	public static double upiupig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.5 * (
				+ f1g2d(0, 0, 1, -1, D11, D12, p11, p12, R, A, B, tau1, tau2)
						- f1g2d(0, 0, 1, +1, D11, D12, p11, p12, R, A, B, tau1, tau2)
		);
	}

	public static double uzuzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.25 * (
				+ f1g2d(+1, -1, 0, 0, D11, D12, p11, p12, R, A, B, tau1, tau2)
						+ f1g2d(-1, +1, 0, 0, D11, D12, p11, p12, R, A, B, tau1, tau2)
						- f1g2d(-1, -1, 0, 0, D11, D12, p11, p12, R, A, B, tau1, tau2)
						- f1g2d(+1, +1, 0, 0, D11, D12, p11, p12, R, A, B, tau1, tau2)
		);
	}

	public static double upiQpizg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.25 * (
				+ f1g2d(0, -1, +1, +1, D11, D22, p11, p22, R, A, B, tau1, tau2)
						- f1g2d(0, -1, +1, -1, D11, D22, p11, p22, R, A, B, tau1, tau2)
						+ f1g2d(0, +1, +1, -1, D11, D22, p11, p22, R, A, B, tau1, tau2)
						- f1g2d(0, +1, +1, +1, D11, D22, p11, p22, R, A, B, tau1, tau2)
		);
	}

	public static double uzQpipig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.25 * (
				- f1g2d(+1, 0, 0, 2, D11, D22, p11, p22, R, A, B, tau1, tau2)
						+ f1g2d(-1, 0, 0, 2, D11, D22, p11, p22, R, A, B, tau1, tau2)
						+ f1g2d(+1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2)
						- f1g2d(-1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2)
		);
	}

	public static double uzQzzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.125 * (
				- f1g2d(+1, -2, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2)
						+ f1g2d(-1, -2, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2)
						- f1g2d(+1, +2, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2)
						+ f1g2d(-1, +2, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2)
						+ 2 * f1g2d(+1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2)
						- 2 * f1g2d(-1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2)
		);
	}

	public static double QpipiQpipig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.125 * (
				+ f1g2d(0, 0, +2, -2, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ f1g2d(0, 0, +2, +2, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- 2 * f1g2d(0, 0, +2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- 2 * f1g2d(0, 0, 0, +2, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ 2 * f1g2d(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
		);
	}

	public static double QxxQyyg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.25 * (
				+ f2g2d(0, 0, 2, 2, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- f2g2d(0, 0, 2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- f2g2d(0, 0, 0, 2, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ f2g2d(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
		);
	}

	public static double QpipiQzzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.125 * (
				+ f1g2d(0, -2, +2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ f1g2d(0, +2, +2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- f1g2d(0, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- f1g2d(0, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- 2 * f1g2d(0, 0, 2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ 2 * f1g2d(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
		);
	}


	public static double QzzQzzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.0625 * (
				+ f1g2d(+2, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ f1g2d(+2, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ f1g2d(-2, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ f1g2d(-2, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- 2 * f1g2d(+2, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- 2 * f1g2d(-2, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- 2 * f1g2d(0, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- 2 * f1g2d(0, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ 4 * f1g2d(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2)
		);
	}

	public static double QpizQpizg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.125 * (
				+ f1g2d(+1, -1, +1, -1, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- f1g2d(+1, -1, +1, +1, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- f1g2d(+1, +1, +1, -1, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ f1g2d(+1, +1, +1, +1, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- f1g2d(-1, -1, +1, -1, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ f1g2d(-1, -1, +1, +1, D21, D22, p21, p22, R, A, B, tau1, tau2)
						+ f1g2d(-1, +1, +1, -1, D21, D22, p21, p22, R, A, B, tau1, tau2)
						- f1g2d(-1, +1, +1, +1, D21, D22, p21, p22, R, A, B, tau1, tau2)
		);
	}

	public static double ssssg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ssppippig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double sspzpzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ppippissg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double pzpzssg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double ppippippippig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) + QpipiQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double pxpxpypyg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) + QxxQyyg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ppippipzpzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) + QpipiQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double pzpzppippig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) + QpipiQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double pzpzpzpzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) + qQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) + QzzQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double spzssg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return -quzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double spzppippig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return -quzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) + uzQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double spzpzpzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return -quzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) + uzQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ssspzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return quzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ppippispzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return quzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) - uzQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double pzpzspzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return quzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) - uzQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double sppisppig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return upiupig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double spzspzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return uzuzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double sppippipzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return upiQpizg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ppipzsppig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return -upiQpizg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double ppipzppipzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return QpizQpizg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double pxpypxpyg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, double[] A, double[] B, int tau1, int tau2) {
		return 0.5 * (ppippippippig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) - pxpxpypyg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2));
	}

	public static double qqpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0;
	}

	public static double g0 (int a2, int b2, double Dn2, double p01, double pn2, double R, double Dn2pd, double pn2pd) {
		return a2 * Dn2pd * (R + a2 * Dn2) + b2 * b2 * Dn2 * Dn2pd + (p01 + pn2) * pn2pd;
	}

	public static double f0pdpartial (int a2, int b2, double Dn2, double p01, double pn2, double R, double Dn2pd, double pn2pd) {

		double f0 = f0(a2, b2, Dn2, p01, pn2, R);

		return  - g0(a2, b2, Dn2, p01, pn2, R, Dn2pd, pn2pd) * Math.pow(f0, 3);
	}

	public static double f0pd (int a2, int b2, double Dn2, double p01, double pn2, double R, int num, double Dn2pd, double pn2pd) {

		switch (num) {
			case 0:
				return 0;
			case 1:
			case 2:
				return f0pdpartial(a2, b2,Dn2, p01, pn2, R, Dn2pd, pn2pd);
		}

		return 0;
	}

	public static double quzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.5 * (
				+ f0pd(+1, 0, D12, p01, p12, R, num, D1pd, p1pd)
						- f0pd(-1, 0, D12, p01, p12, R, num, D1pd, p1pd)
		);
	}

	public static double qQpipipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.5 * (
				+ f0pd(0, 2, D22, p01, p22, R, num, D2pd, p2pd)
						- f0pd(0, 0, D22, p01, p22, R, num, D2pd, p2pd)
		);

	}

	public static double qQzzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.25 * (
				+ f0pd(+2, 0, D22, p01, p22, R, num, D2pd, p2pd)
						+ f0pd(-2, 0, D22, p01, p22, R, num, D2pd, p2pd)
						- 2 * f0pd(0, 0, D22, p01, p22, R, num, D2pd, p2pd)
		);

	}

	public static double g1 (int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R, double Dxpd, double pxpd, int type) {
		switch (type) {
			case 0:
				return Dxpd * (a1 * (R + a1 * Dn1 + a2 * Dm2) + b1 * (b1 * Dn1 + b2 * Dm2)) + (pn1 + pm2) * pxpd;
			case 1:
				return Dxpd * (a2 * (R + a1 * Dn1 + a2 * Dm2) + b2 * (b1 * Dn1 + b2 * Dm2)) + (pn1 + pm2) * pxpd;
		}
		return 0;
	}

	public static double f1pd (int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R, int type, double Dxpd, double pxpd) {

		switch (type) {
			case 0:
			case 1:
				return - g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, type) * Math.pow(f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
			case 2:
				double f1 = f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);

				return - g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, 0) * Math.pow(f1, 3)
						- g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, 1) * Math.pow(f1, 3);
		}

		return 0;

	}

	public static double g2(int a1, int a2, int b1, int b2, double D21, double D22, double p21, double p22, double R, double D2pd, double p2pd, int type) {
		switch (type) {
			case 0:
				return D2pd * (a1 * (R + a1 * D21 + a2 * D22) + b1 * b1 * D21) + (p21 + p22) * p2pd;
			case 1:
				return D2pd * (a2 * (R + a1 * D21 + a2 * D22) + b2 * b2 * D22) + (p21 + p22) * p2pd;
		}
		return 0;
	}


	public static double f2pd (int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R, int type, double Dxpd, double pxpd) {

		switch (type) {
			case 0:
			case 1:
				return - g2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, type) * Math.pow(f2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
			case 2:
				double f2 = f2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);
				return - g2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, 0) * Math.pow(f2, 3)
						- g2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, 1) * Math.pow(f2, 3);
		}

		return 0;

	}

	public static double upiupipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.5 * (
				+ f1pd(0, 0, 1, -1, D11, D12, p11, p12, R, num, D1pd, p1pd)
						- f1pd(0, 0, 1, +1, D11, D12, p11, p12, R, num, D1pd, p1pd)
		);

	}

	public static double uzuzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.25 * (
				+ f1pd(+1, -1, 0, 0, D11, D12, p11, p12, R, num, D1pd, p1pd)
						+ f1pd(-1, +1, 0, 0, D11, D12, p11, p12, R, num, D1pd, p1pd)
						- f1pd(-1, -1, 0, 0, D11, D12, p11, p12, R, num, D1pd, p1pd)
						- f1pd(+1, +1, 0, 0, D11, D12, p11, p12, R, num, D1pd, p1pd)
		);

	}

	public static double varf1pd (int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R, int type, double Dn1pd, double pn1pd, double Dm2pd, double pm2pd) {

		switch (type) {
			case 0:
				return - g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dn1pd, pn1pd, type) * Math.pow(f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
			case 1:
				return - g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dm2pd, pm2pd, type) * Math.pow(f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
			case 2:
				double f1 = f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);
				return - g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dn1pd, pn1pd, 0) * Math.pow(f1, 3)
						- g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dm2pd, pm2pd, 1) * Math.pow(f1, 3);
		}

		return 0;

	}

	public static double upiQpizpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.25 * (
				+ varf1pd(0, -1, +1, +1, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						- varf1pd(0, -1, +1, -1, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						+ varf1pd(0, +1, +1, -1, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						- varf1pd(0, +1, +1, +1, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
		);

	}

	public static double uzQpipipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.25 * (
				- varf1pd(+1, 0, 0, 2, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						+ varf1pd(-1, 0, 0, 2, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						+ varf1pd(+1, 0, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						- varf1pd(-1, 0, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
		);

	}

	public static double uzQzzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.125 * (
				- varf1pd(+1, -2, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						+ varf1pd(-1, -2, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						- varf1pd(+1, +2, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						+ varf1pd(-1, +2, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						+ 2 * varf1pd(+1, 0, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
						- 2 * varf1pd(-1, 0, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd)
		);

	}

	public static double QpipiQpipipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.125 * (
				+ f1pd(0, 0, +2, -2, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ f1pd(0, 0, +2, +2, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- 2 * f1pd(0, 0, +2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- 2 * f1pd(0, 0, 0, +2, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ 2 * f1pd(0, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
		);
	}

	public static double QxxQyypd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.25 * (
				+ f2pd(0, 0, 2, 2, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- f2pd(0, 0, 2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- f2pd(0, 0, 0, 2, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ f2pd(0, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
		);

	}

	public static double QpipiQzzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.125 * (
				+ f1pd(0, -2, +2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ f1pd(0, +2, +2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- f1pd(0, -2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- f1pd(0, +2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- 2 * f1pd(0, 0, 2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ 2 * f1pd(0, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
		);

	}

	public static double QzzQzzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.0625 * (
				+ f1pd(+2, -2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ f1pd(+2, +2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ f1pd(-2, -2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ f1pd(-2, +2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- 2 * f1pd(+2, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- 2 * f1pd(-2, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- 2 * f1pd(0, +2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- 2 * f1pd(0, -2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ 4 * f1pd(0, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd)
		);
	}

	public static double QpizQpizpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double R, int num, double D1pd, double D2pd, double p1pd, double p2pd) {
		return 0.125 * (
				+ f1pd(+1, -1, +1, -1, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- f1pd(+1, -1, +1, +1, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- f1pd(+1, +1, +1, -1, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ f1pd(+1, +1, +1, +1, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- f1pd(-1, -1, +1, -1, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ f1pd(-1, -1, +1, +1, D21, D22, p21, p22, R, num, D2pd, p2pd)
						+ f1pd(-1, +1, +1, -1, D21, D22, p21, p22, R, num, D2pd, p2pd)
						- f1pd(-1, +1, +1, +1, D21, D22, p21, p22, R, num, D2pd, p2pd)
		);
	}

	public static double sssspd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, int num, double D1deriv,
								   double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ssppippipd(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double sspzpzpd(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippisspd(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzsspd(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippippippipd(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, int num,
										   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv) +
				QpipiQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pxpxpypypd(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv) +
				QxxQyypd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippipzpzpd(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv) +
				QpipiQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double pzpzppippipd(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				QpipiQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzpzpzpd(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				QzzQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzsspd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return -quzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzppippipd(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return -quzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				uzQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzpzpzpd(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, int num,
									  double D1deriv,
									  double D2deriv, double p1deriv, double p2deriv) {
		return -quzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				uzQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ssspzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return quzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippispzpd(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return quzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) -
				uzQpipipd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzspzpd(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, int num,
									  double D1deriv,
									  double D2deriv, double p1deriv, double p2deriv) {
		return quzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) -
				uzQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double sppisppipd(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return upiupipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzspzpd(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return uzuzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double sppippipzpd(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return upiQpizpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppipzsppipd(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return -upiQpizpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
				p2deriv);
	}

	public static double ppipzppipzpd(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return QpizQpizpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double pxpypxpypd(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return 0.5 *
				(ppippippippipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv) -
						pxpxpypypd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv,
								p1deriv, p2deriv));
	}

	public static int f(int num) {
		if (num == 0 || num == 1) return 1 - num;
		else if (num == 2) return num;
		return -1;
	}

}
