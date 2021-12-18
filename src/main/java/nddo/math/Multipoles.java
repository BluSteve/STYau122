package nddo.math;

public class Multipoles {
	public static double qq(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							double p22, double D12, double D22, double R) {
		return Math.pow(R * R + (p01 + p02) * (p01 + p02), -0.5);
	}

	public static double f0(int a2, int b2, double Dn2, double p01, double pn2, double R) {
		return Math.pow((R + a2 * Dn2) * (R + a2 * Dn2) + b2 * b2 * Dn2 * Dn2 + (p01 + pn2) * (p01 + pn2), -0.5);
	}

	public static double quz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							 double p22, double D12, double D22, double R) {
		return 0.5 * (+f0(+1, 0, D12, p01, p12, R) - f0(-1, 0, D12, p01, p12, R));
	}

	public static double qQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		return 0.5 * (+f0(0, 2, D22, p01, p22, R) - f0(0, 0, D22, p01, p22, R));
	}

	public static double qQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							  double p22, double D12, double D22, double R) {
		return 0.25 * (+f0(+2, 0, D22, p01, p22, R) + f0(-2, 0, D22, p01, p22, R) - 2 * f0(0, 0, D22, p01, p22, R));
	}

	public static double f1(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R) {
		return Math.pow(
				(R + a1 * Dn1 + a2 * Dm2) * (R + a1 * Dn1 + a2 * Dm2) + (b1 * Dn1 + b2 * Dm2) * (b1 * Dn1 + b2 * Dm2) +
						(pn1 + pm2) * (pn1 + pm2), -0.5);
	}

	public static double f2(int a1, int a2, int b1, int b2, double D21, double D22, double p21, double p22, double R) {
		return Math.pow((R + a1 * D21 + a2 * D22) * (R + a1 * D21 + a2 * D22) + (b1 * D21) * (b1 * D21) +
				(b2 * D22) * (b2 * D22) + (p21 + p22) * (p21 + p22), -0.5);
	}

	public static double upiupi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		return 0.5 * (+f1(0, 0, 1, -1, D11, D12, p11, p12, R) - f1(0, 0, 1, +1, D11, D12, p11, p12, R));
	}

	public static double uzuz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							  double p22, double D12, double D22, double R) {
		return 0.25 * (+f1(+1, -1, 0, 0, D11, D12, p11, p12, R) + f1(-1, +1, 0, 0, D11, D12, p11, p12, R) -
				f1(-1, -1, 0, 0, D11, D12, p11, p12, R) - f1(+1, +1, 0, 0, D11, D12, p11, p12, R));
	}

	public static double upiQpiz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R) {
		return 0.25 * (+f1(0, -1, +1, +1, D11, D22, p11, p22, R) - f1(0, -1, +1, -1, D11, D22, p11, p22, R) +
				f1(0, +1, +1, -1, D11, D22, p11, p22, R) - f1(0, +1, +1, +1, D11, D22, p11, p22, R));
	}

	public static double uzQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R) {
		return 0.25 * (-f1(+1, 0, 0, 2, D11, D22, p11, p22, R) + f1(-1, 0, 0, 2, D11, D22, p11, p22, R) +
				f1(+1, 0, 0, 0, D11, D22, p11, p22, R) - f1(-1, 0, 0, 0, D11, D22, p11, p22, R));
	}

	public static double uzQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							   double p22, double D12, double D22, double R) {
		return 0.125 * (-f1(+1, -2, 0, 0, D11, D22, p11, p22, R) + f1(-1, -2, 0, 0, D11, D22, p11, p22, R) -
				f1(+1, +2, 0, 0, D11, D22, p11, p22, R) + f1(-1, +2, 0, 0, D11, D22, p11, p22, R) +
				2 * f1(+1, 0, 0, 0, D11, D22, p11, p22, R) - 2 * f1(-1, 0, 0, 0, D11, D22, p11, p22, R));
	}

	public static double QpipiQpipi(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R) {
		return 0.125 * (+f1(0, 0, +2, -2, D21, D22, p21, p22, R) + f1(0, 0, +2, +2, D21, D22, p21, p22, R) -
				2 * f1(0, 0, +2, 0, D21, D22, p21, p22, R) - 2 * f1(0, 0, 0, +2, D21, D22, p21, p22, R) +
				2 * f1(0, 0, 0, 0, D21, D22, p21, p22, R));
	}

	public static double QxxQyy(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		return 0.25 * (+f2(0, 0, 2, 2, D21, D22, p21, p22, R) - f2(0, 0, 2, 0, D21, D22, p21, p22, R) -
				f2(0, 0, 0, 2, D21, D22, p21, p22, R) + f2(0, 0, 0, 0, D21, D22, p21, p22, R));
	}

	public static double QpipiQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		return 0.125 * (+f1(0, -2, +2, 0, D21, D22, p21, p22, R) + f1(0, +2, +2, 0, D21, D22, p21, p22, R) -
				f1(0, -2, 0, 0, D21, D22, p21, p22, R) - f1(0, +2, 0, 0, D21, D22, p21, p22, R) -
				2 * f1(0, 0, 2, 0, D21, D22, p21, p22, R) + 2 * f1(0, 0, 0, 0, D21, D22, p21, p22, R));
	}

	public static double QzzQzz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R) {
		return 0.0625 * (+f1(+2, -2, 0, 0, D21, D22, p21, p22, R) + f1(+2, +2, 0, 0, D21, D22, p21, p22, R) +
				f1(-2, -2, 0, 0, D21, D22, p21, p22, R) + f1(-2, +2, 0, 0, D21, D22, p21, p22, R) -
				2 * f1(+2, 0, 0, 0, D21, D22, p21, p22, R) - 2 * f1(-2, 0, 0, 0, D21, D22, p21, p22, R) -
				2 * f1(0, +2, 0, 0, D21, D22, p21, p22, R) - 2 * f1(0, -2, 0, 0, D21, D22, p21, p22, R) +
				4 * f1(0, 0, 0, 0, D21, D22, p21, p22, R));
	}

	public static double QpizQpiz(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R) {
		return 0.125 * (+f1(+1, -1, +1, -1, D21, D22, p21, p22, R) - f1(+1, -1, +1, +1, D21, D22, p21, p22, R) -
				f1(+1, +1, +1, -1, D21, D22, p21, p22, R) + f1(+1, +1, +1, +1, D21, D22, p21, p22, R) -
				f1(-1, -1, +1, -1, D21, D22, p21, p22, R) + f1(-1, -1, +1, +1, D21, D22, p21, p22, R) +
				f1(-1, +1, +1, -1, D21, D22, p21, p22, R) - f1(-1, +1, +1, +1, D21, D22, p21, p22, R));
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

	public static double qqgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							  double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return (B[tau] - A[tau]) * Math.pow(R * R + (p01 + p02) * (p01 + p02), -1.5);
	}

	public static double f0gd(int a2, int b2, double Dn2, double p01, double pn2, double R, double[] A, double[] B,
							  int tau) {
		return (B[tau] - A[tau]) / R * (R + a2 * Dn2) * Math.pow(f0(a2, b2, Dn2, p01, pn2, R), 3);
	}

	public static double quzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							   double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.5 * (+f0gd(+1, 0, D12, p01, p12, R, A, B, tau) - f0gd(-1, 0, D12, p01, p12, R, A, B, tau));
	}

	public static double qQpipigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.5 * (+f0gd(0, 2, D22, p01, p22, R, A, B, tau) - f0gd(0, 0, D22, p01, p22, R, A, B, tau));
	}

	public static double qQzzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (+f0gd(+2, 0, D22, p01, p22, R, A, B, tau) + f0gd(-2, 0, D22, p01, p22, R, A, B, tau) -
				2 * f0gd(0, 0, D22, p01, p22, R, A, B, tau));
	}

	public static double f1gd(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R,
							  double[] A, double[] B, int tau) {
		return (B[tau] - A[tau]) / R * (R + a1 * Dn1 + a2 * Dm2) *
				Math.pow(f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
	}

	public static double f2gd(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R,
							  double[] A, double[] B, int tau) {
		return (B[tau] - A[tau]) / R * (R + a1 * Dn1 + a2 * Dm2) *
				Math.pow(f2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
	}

	public static double upiupigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.5 * (+f1gd(0, 0, 1, -1, D11, D12, p11, p12, R, A, B, tau) -
				f1gd(0, 0, 1, +1, D11, D12, p11, p12, R, A, B, tau));
	}

	public static double uzuzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (+f1gd(+1, -1, 0, 0, D11, D12, p11, p12, R, A, B, tau) +
				f1gd(-1, +1, 0, 0, D11, D12, p11, p12, R, A, B, tau) -
				f1gd(-1, -1, 0, 0, D11, D12, p11, p12, R, A, B, tau) -
				f1gd(+1, +1, 0, 0, D11, D12, p11, p12, R, A, B, tau));
	}

	public static double upiQpizgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (+f1gd(0, -1, +1, +1, D11, D22, p11, p22, R, A, B, tau) -
				f1gd(0, -1, +1, -1, D11, D22, p11, p22, R, A, B, tau) +
				f1gd(0, +1, +1, -1, D11, D22, p11, p22, R, A, B, tau) -
				f1gd(0, +1, +1, +1, D11, D22, p11, p22, R, A, B, tau));
	}

	public static double uzQpipigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (-f1gd(+1, 0, 0, 2, D11, D22, p11, p22, R, A, B, tau) +
				f1gd(-1, 0, 0, 2, D11, D22, p11, p22, R, A, B, tau) +
				f1gd(+1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau) -
				f1gd(-1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau));
	}

	public static double uzQzzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.125 * (-f1gd(+1, -2, 0, 0, D11, D22, p11, p22, R, A, B, tau) +
				f1gd(-1, -2, 0, 0, D11, D22, p11, p22, R, A, B, tau) -
				f1gd(+1, +2, 0, 0, D11, D22, p11, p22, R, A, B, tau) +
				f1gd(-1, +2, 0, 0, D11, D22, p11, p22, R, A, B, tau) +
				2 * f1gd(+1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau) -
				2 * f1gd(-1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau));
	}

	public static double QpipiQpipigd(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, double[] A, double[] B,
									  int tau) {
		return 0.125 * (+f1gd(0, 0, +2, -2, D21, D22, p21, p22, R, A, B, tau) +
				f1gd(0, 0, +2, +2, D21, D22, p21, p22, R, A, B, tau) -
				2 * f1gd(0, 0, +2, 0, D21, D22, p21, p22, R, A, B, tau) -
				2 * f1gd(0, 0, 0, +2, D21, D22, p21, p22, R, A, B, tau) +
				2 * f1gd(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau));
	}

	public static double QxxQyygd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.25 * (+f2gd(0, 0, 2, 2, D21, D22, p21, p22, R, A, B, tau) -
				f2gd(0, 0, 2, 0, D21, D22, p21, p22, R, A, B, tau) -
				f2gd(0, 0, 0, 2, D21, D22, p21, p22, R, A, B, tau) +
				f2gd(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau));
	}

	public static double QpipiQzzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.125 * (+f1gd(0, -2, +2, 0, D21, D22, p21, p22, R, A, B, tau) +
				f1gd(0, +2, +2, 0, D21, D22, p21, p22, R, A, B, tau) -
				f1gd(0, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau) -
				f1gd(0, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau) -
				2 * f1gd(0, 0, 2, 0, D21, D22, p21, p22, R, A, B, tau) +
				2 * f1gd(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau));
	}

	public static double QzzQzzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.0625 * (+f1gd(+2, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau) +
				f1gd(+2, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau) +
				f1gd(-2, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau) +
				f1gd(-2, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau) -
				2 * f1gd(+2, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau) -
				2 * f1gd(-2, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau) -
				2 * f1gd(0, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau) -
				2 * f1gd(0, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau) +
				4 * f1gd(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau));
	}

	public static double QpizQpizgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.125 * (+f1gd(+1, -1, +1, -1, D21, D22, p21, p22, R, A, B, tau) -
				f1gd(+1, -1, +1, +1, D21, D22, p21, p22, R, A, B, tau) -
				f1gd(+1, +1, +1, -1, D21, D22, p21, p22, R, A, B, tau) +
				f1gd(+1, +1, +1, +1, D21, D22, p21, p22, R, A, B, tau) -
				f1gd(-1, -1, +1, -1, D21, D22, p21, p22, R, A, B, tau) +
				f1gd(-1, -1, +1, +1, D21, D22, p21, p22, R, A, B, tau) +
				f1gd(-1, +1, +1, -1, D21, D22, p21, p22, R, A, B, tau) -
				f1gd(-1, +1, +1, +1, D21, D22, p21, p22, R, A, B, tau));
	}

	public static double ssssgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ssppippigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double sspzpzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ppippissgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double pzpzssgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double ppippippippigd(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, double[] A,
										double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) +
				QpipiQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double pxpxpypygd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) +
				QxxQyygd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ppippipzpzgd(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, double[] A, double[] B,
									  int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) +
				QpipiQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double pzpzppippigd(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, double[] A, double[] B,
									  int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) +
				QpipiQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double pzpzpzpzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return qqgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) +
				qQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) +
				QzzQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double spzssgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return -quzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double spzppippigd(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return -quzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) +
				uzQpipigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double spzpzpzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return -quzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau) +
				uzQzzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ssspzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return quzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ppippispzgd(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return quzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) -
				uzQpipigd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double pzpzspzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return quzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) -
				uzQzzgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double sppisppigd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return upiupigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double spzspzgd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return uzuzgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double sppippipzgd(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return upiQpizgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double ppipzsppigd(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return -upiQpizgd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau);
	}

	public static double ppipzppipzgd(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, double[] A, double[] B,
									  int tau) {
		return QpizQpizgd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau);
	}

	public static double pxpypxpygd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau) {
		return 0.5 * (ppippippippigd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau) -
				pxpxpypygd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau));
	}

	public static double qqg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							   double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
							   int tau2) {
		double sum = 0;
		if (tau1 == tau2) sum = -Math.pow(R * R + (p01 + p02) * (p01 + p02), -1.5);
		return sum + 3 * (B[tau1] - A[tau1]) * (B[tau2] - A[tau2]) * Math.pow(R * R + (p01 + p02) * (p01 + p02), -2.5);
	}

	public static double f0g2d(int a2, int b2, double Dn2, double p01, double pn2, double R, double[] A, double[] B,
							   int tau1, int tau2) {
		double f0 = f0(a2, b2, Dn2, p01, pn2, R);
		double num = (B[tau1] - A[tau1]) * (B[tau2] - A[tau2]) / (R * R) * Math.pow(f0, 3) *
				(3 * Math.pow(R + a2 * Dn2, 2) * f0 * f0 + a2 * Dn2 / R);
		if (tau1 == tau2) num += -(R + a2 * Dn2) * Math.pow(f0, 3) / R;
		return num;
	}

	public static double quzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								int tau2) {
		return 0.5 *
				(+f0g2d(+1, 0, D12, p01, p12, R, A, B, tau1, tau2) - f0g2d(-1, 0, D12, p01, p12, R, A, B, tau1, tau2));
	}

	public static double qQpipig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								   int tau2) {
		return 0.5 *
				(+f0g2d(0, 2, D22, p01, p22, R, A, B, tau1, tau2) - f0g2d(0, 0, D22, p01, p22, R, A, B, tau1, tau2));
	}

	public static double qQzzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								 int tau2) {
		return 0.25 *
				(+f0g2d(+2, 0, D22, p01, p22, R, A, B, tau1, tau2) + f0g2d(-2, 0, D22, p01, p22, R, A, B, tau1, tau2) -
						2 * f0g2d(0, 0, D22, p01, p22, R, A, B, tau1, tau2));
	}

	public static double f1g2d(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2,
							   double R,
							   double[] A, double[] B, int tau1, int tau2) {
		double f1 = f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);
		double num = (B[tau1] - A[tau1]) * (B[tau2] - A[tau2]) / (R * R) * Math.pow(f1, 3) *
				(3 * Math.pow(R + a1 * Dn1 + a2 * Dm2, 2) * f1 * f1 + (a1 * Dn1 + a2 * Dm2) / R);
		if (tau1 == tau2) num += -(R + a1 * Dn1 + a2 * Dm2) * Math.pow(f1, 3) / R;
		return num;
	}

	public static double f2g2d(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2,
							   double R,
							   double[] A, double[] B, int tau1, int tau2) {
		double f2 = f2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);
		double num = (B[tau1] - A[tau1]) * (B[tau2] - A[tau2]) / (R * R) * Math.pow(f2, 3) *
				(3 * Math.pow(R + a1 * Dn1 + a2 * Dm2, 2) * f2 * f2 + (a1 * Dn1 + a2 * Dm2) / R);
		if (tau1 == tau2) num += -(R + a1 * Dn1 + a2 * Dm2) * Math.pow(f2, 3) / R;
		return num;
	}

	public static double upiupig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								   int tau2) {
		return 0.5 * (+f1g2d(0, 0, 1, -1, D11, D12, p11, p12, R, A, B, tau1, tau2) -
				f1g2d(0, 0, 1, +1, D11, D12, p11, p12, R, A, B, tau1, tau2));
	}

	public static double uzuzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								 int tau2) {
		return 0.25 * (+f1g2d(+1, -1, 0, 0, D11, D12, p11, p12, R, A, B, tau1, tau2) +
				f1g2d(-1, +1, 0, 0, D11, D12, p11, p12, R, A, B, tau1, tau2) -
				f1g2d(-1, -1, 0, 0, D11, D12, p11, p12, R, A, B, tau1, tau2) -
				f1g2d(+1, +1, 0, 0, D11, D12, p11, p12, R, A, B, tau1, tau2));
	}

	public static double upiQpizg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									int tau2) {
		return 0.25 * (+f1g2d(0, -1, +1, +1, D11, D22, p11, p22, R, A, B, tau1, tau2) -
				f1g2d(0, -1, +1, -1, D11, D22, p11, p22, R, A, B, tau1, tau2) +
				f1g2d(0, +1, +1, -1, D11, D22, p11, p22, R, A, B, tau1, tau2) -
				f1g2d(0, +1, +1, +1, D11, D22, p11, p22, R, A, B, tau1, tau2));
	}

	public static double uzQpipig2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									int tau2) {
		return 0.25 * (-f1g2d(+1, 0, 0, 2, D11, D22, p11, p22, R, A, B, tau1, tau2) +
				f1g2d(-1, 0, 0, 2, D11, D22, p11, p22, R, A, B, tau1, tau2) +
				f1g2d(+1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2) -
				f1g2d(-1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2));
	}

	public static double uzQzzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								  int tau2) {
		return 0.125 * (-f1g2d(+1, -2, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2) +
				f1g2d(-1, -2, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2) -
				f1g2d(+1, +2, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2) +
				f1g2d(-1, +2, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2) +
				2 * f1g2d(+1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2) -
				2 * f1g2d(-1, 0, 0, 0, D11, D22, p11, p22, R, A, B, tau1, tau2));
	}

	public static double QpipiQpipig2d(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, double[] A,
									   double[] B,
									   int tau1, int tau2) {
		return 0.125 * (+f1g2d(0, 0, +2, -2, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				f1g2d(0, 0, +2, +2, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				2 * f1g2d(0, 0, +2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				2 * f1g2d(0, 0, 0, +2, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				2 * f1g2d(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2));
	}

	public static double QxxQyyg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								   int tau2) {
		return 0.25 * (+f2g2d(0, 0, 2, 2, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				f2g2d(0, 0, 2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				f2g2d(0, 0, 0, 2, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				f2g2d(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2));
	}

	public static double QpipiQzzg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									 int tau2) {
		return 0.125 * (+f1g2d(0, -2, +2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				f1g2d(0, +2, +2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				f1g2d(0, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				f1g2d(0, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				2 * f1g2d(0, 0, 2, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				2 * f1g2d(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2));
	}

	public static double QzzQzzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								   int tau2) {
		return 0.0625 * (+f1g2d(+2, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				f1g2d(+2, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				f1g2d(-2, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				f1g2d(-2, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				2 * f1g2d(+2, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				2 * f1g2d(-2, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				2 * f1g2d(0, +2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				2 * f1g2d(0, -2, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				4 * f1g2d(0, 0, 0, 0, D21, D22, p21, p22, R, A, B, tau1, tau2));
	}

	public static double QpizQpizg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									 int tau2) {
		return 0.125 * (+f1g2d(+1, -1, +1, -1, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				f1g2d(+1, -1, +1, +1, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				f1g2d(+1, +1, +1, -1, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				f1g2d(+1, +1, +1, +1, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				f1g2d(-1, -1, +1, -1, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				f1g2d(-1, -1, +1, +1, D21, D22, p21, p22, R, A, B, tau1, tau2) +
				f1g2d(-1, +1, +1, -1, D21, D22, p21, p22, R, A, B, tau1, tau2) -
				f1g2d(-1, +1, +1, +1, D21, D22, p21, p22, R, A, B, tau1, tau2));
	}

	public static double ssssg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								 int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ssppippig2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									 int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double sspzpzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								   int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ppippissg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									 int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double pzpzssg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								   int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double ppippippippig2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, double[] A,
										 double[] B, int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) +
				QpipiQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double pxpxpypyg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									 int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) +
				QxxQyyg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ppippipzpzg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, double[] A,
									   double[] B,
									   int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) +
				QpipiQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double pzpzppippig2d(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, double[] A,
									   double[] B,
									   int tau1, int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) +
				QpipiQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double pzpzpzpzg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									 int tau2) {
		return qqg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) +
				qQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) +
				QzzQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double spzssg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								  int tau2) {
		return -quzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double spzppippig2d(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, double[] A, double[] B,
									  int tau1, int tau2) {
		return -quzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) +
				uzQpipig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double spzpzpzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									int tau2) {
		return -quzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2) +
				uzQzzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ssspzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								  int tau2) {
		return quzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ppippispzg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, double[] A, double[] B,
									  int tau1, int tau2) {
		return quzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) -
				uzQpipig2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double pzpzspzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									int tau2) {
		return quzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) -
				uzQzzg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double sppisppig2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									 int tau2) {
		return upiupig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double spzspzg2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
								   int tau2) {
		return uzuzg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double sppippipzg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, double[] A, double[] B,
									  int tau1, int tau2) {
		return upiQpizg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double ppipzsppig2d(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, double[] A, double[] B,
									  int tau1, int tau2) {
		return -upiQpizg2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, A, B, tau1, tau2);
	}

	public static double ppipzppipzg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, double[] A,
									   double[] B,
									   int tau1, int tau2) {
		return QpizQpizg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2);
	}

	public static double pxpypxpyg2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double[] A, double[] B, int tau1,
									 int tau2) {
		return 0.5 * (ppippippippig2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2) -
				pxpxpypyg2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, A, B, tau1, tau2));
	}

	public static double qqpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							  double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
							  double p1pd, double p2pd) {
		return 0;
	}

	public static double g0(int a2, int b2, double Dn2, double p01, double pn2, double R, double Dn2pd, double pn2pd) {
		return a2 * Dn2pd * (R + a2 * Dn2) + b2 * b2 * Dn2 * Dn2pd + (p01 + pn2) * pn2pd;
	}

	public static double f0pdpartial(int a2, int b2, double Dn2, double p01, double pn2, double R, double Dn2pd,
									 double pn2pd) {
		double f0 = f0(a2, b2, Dn2, p01, pn2, R);
		return -g0(a2, b2, Dn2, p01, pn2, R, Dn2pd, pn2pd) * Math.pow(f0, 3);
	}

	public static double f0pd(int a2, int b2, double Dn2, double p01, double pn2, double R, int num, double Dn2pd,
							  double pn2pd) {
		switch (num) {
			case 0:
				return 0;
			case 1:
			case 2:
				return f0pdpartial(a2, b2, Dn2, p01, pn2, R, Dn2pd, pn2pd);
		}
		return 0;
	}

	public static double quzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
							   double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
							   double p1pd, double p2pd) {
		return 0.5 * (+f0pd(+1, 0, D12, p01, p12, R, num, D1pd, p1pd) - f0pd(-1, 0, D12, p01, p12, R, num, D1pd,
				p1pd));
	}

	public static double qQpipipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
								  double p1pd, double p2pd) {
		return 0.5 * (+f0pd(0, 2, D22, p01, p22, R, num, D2pd, p2pd) - f0pd(0, 0, D22, p01, p22, R, num, D2pd, p2pd));
	}

	public static double qQzzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
								double p1pd, double p2pd) {
		return 0.25 *
				(+f0pd(+2, 0, D22, p01, p22, R, num, D2pd, p2pd) + f0pd(-2, 0, D22, p01, p22, R, num, D2pd, p2pd) -
						2 * f0pd(0, 0, D22, p01, p22, R, num, D2pd, p2pd));
	}

	public static double g1(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R,
							double Dxpd, double pxpd, int type) {
		switch (type) {
			case 0:
				return Dxpd * (a1 * (R + a1 * Dn1 + a2 * Dm2) + b1 * (b1 * Dn1 + b2 * Dm2)) + (pn1 + pm2) * pxpd;
			case 1:
				return Dxpd * (a2 * (R + a1 * Dn1 + a2 * Dm2) + b2 * (b1 * Dn1 + b2 * Dm2)) + (pn1 + pm2) * pxpd;
		}
		return 0;
	}

	public static double f1pd(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R,
							  int type, double Dxpd, double pxpd) {
		switch (type) {
			case 0:
			case 1:
				return -g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, type) *
						Math.pow(f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
			case 2:
				double f1 = f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);
				return -g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, 0) * Math.pow(f1, 3) -
						g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, 1) * Math.pow(f1, 3);
		}
		return 0;
	}

	public static double g2(int a1, int a2, int b1, int b2, double D21, double D22, double p21, double p22, double R,
							double D2pd, double p2pd, int type) {
		switch (type) {
			case 0:
				return D2pd * (a1 * (R + a1 * D21 + a2 * D22) + b1 * b1 * D21) + (p21 + p22) * p2pd;
			case 1:
				return D2pd * (a2 * (R + a1 * D21 + a2 * D22) + b2 * b2 * D22) + (p21 + p22) * p2pd;
		}
		return 0;
	}

	public static double f2pd(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R,
							  int type, double Dxpd, double pxpd) {
		switch (type) {
			case 0:
			case 1:
				return -g2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, type) *
						Math.pow(f2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
			case 2:
				double f2 = f2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);
				return -g2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, 0) * Math.pow(f2, 3) -
						g2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxpd, pxpd, 1) * Math.pow(f2, 3);
		}
		return 0;
	}

	public static double upiupipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
								  double p1pd, double p2pd) {
		return 0.5 * (+f1pd(0, 0, 1, -1, D11, D12, p11, p12, R, num, D1pd, p1pd) -
				f1pd(0, 0, 1, +1, D11, D12, p11, p12, R, num, D1pd, p1pd));
	}

	public static double uzuzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
								double p1pd, double p2pd) {
		return 0.25 * (+f1pd(+1, -1, 0, 0, D11, D12, p11, p12, R, num, D1pd, p1pd) +
				f1pd(-1, +1, 0, 0, D11, D12, p11, p12, R, num, D1pd, p1pd) -
				f1pd(-1, -1, 0, 0, D11, D12, p11, p12, R, num, D1pd, p1pd) -
				f1pd(+1, +1, 0, 0, D11, D12, p11, p12, R, num, D1pd, p1pd));
	}

	public static double varf1pd(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2,
								 double R, int type, double Dn1pd, double pn1pd, double Dm2pd, double pm2pd) {
		switch (type) {
			case 0:
				return -g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dn1pd, pn1pd, type) *
						Math.pow(f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
			case 1:
				return -g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dm2pd, pm2pd, type) *
						Math.pow(f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R), 3);
			case 2:
				double f1 = f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);
				return -g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dn1pd, pn1pd, 0) * Math.pow(f1, 3) -
						g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dm2pd, pm2pd, 1) * Math.pow(f1, 3);
		}
		return 0;
	}

	public static double upiQpizpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
								   double p1pd, double p2pd) {
		return 0.25 * (+varf1pd(0, -1, +1, +1, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) -
				varf1pd(0, -1, +1, -1, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) +
				varf1pd(0, +1, +1, -1, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) -
				varf1pd(0, +1, +1, +1, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd));
	}

	public static double uzQpipipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
								   double p1pd, double p2pd) {
		return 0.25 * (-varf1pd(+1, 0, 0, 2, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) +
				varf1pd(-1, 0, 0, 2, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) +
				varf1pd(+1, 0, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) -
				varf1pd(-1, 0, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd));
	}

	public static double uzQzzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
								 double p1pd, double p2pd) {
		return 0.125 * (-varf1pd(+1, -2, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) +
				varf1pd(-1, -2, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) -
				varf1pd(+1, +2, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) +
				varf1pd(-1, +2, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) +
				2 * varf1pd(+1, 0, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd) -
				2 * varf1pd(-1, 0, 0, 0, D11, D22, p11, p22, R, num, D1pd, p1pd, D2pd, p2pd));
	}

	public static double QpipiQpipipd(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, int num, double D1pd,
									  double D2pd, double p1pd, double p2pd) {
		return 0.125 * (+f1pd(0, 0, +2, -2, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				f1pd(0, 0, +2, +2, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				2 * f1pd(0, 0, +2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				2 * f1pd(0, 0, 0, +2, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				2 * f1pd(0, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd));
	}

	public static double QxxQyypd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
								  double p1pd, double p2pd) {
		return 0.25 * (+f2pd(0, 0, 2, 2, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				f2pd(0, 0, 2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				f2pd(0, 0, 0, 2, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				f2pd(0, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd));
	}

	public static double QpipiQzzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
									double p1pd, double p2pd) {
		return 0.125 * (+f1pd(0, -2, +2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				f1pd(0, +2, +2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				f1pd(0, -2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				f1pd(0, +2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				2 * f1pd(0, 0, 2, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				2 * f1pd(0, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd));
	}

	public static double QzzQzzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
								  double p1pd, double p2pd) {
		return 0.0625 * (+f1pd(+2, -2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				f1pd(+2, +2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				f1pd(-2, -2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				f1pd(-2, +2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				2 * f1pd(+2, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				2 * f1pd(-2, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				2 * f1pd(0, +2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				2 * f1pd(0, -2, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				4 * f1pd(0, 0, 0, 0, D21, D22, p21, p22, R, num, D2pd, p2pd));
	}

	public static double QpizQpizpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1pd, double D2pd,
									double p1pd, double p2pd) {
		return 0.125 * (+f1pd(+1, -1, +1, -1, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				f1pd(+1, -1, +1, +1, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				f1pd(+1, +1, +1, -1, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				f1pd(+1, +1, +1, +1, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				f1pd(-1, -1, +1, -1, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				f1pd(-1, -1, +1, +1, D21, D22, p21, p22, R, num, D2pd, p2pd) +
				f1pd(-1, +1, +1, -1, D21, D22, p21, p22, R, num, D2pd, p2pd) -
				f1pd(-1, +1, +1, +1, D21, D22, p21, p22, R, num, D2pd, p2pd));
	}

	public static double sssspd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ssppippipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double sspzpzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, int num, double D1deriv,
								  double D2deriv,
								  double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippisspd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzsspd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, int num, double D1deriv,
								  double D2deriv,
								  double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double ppippippippipd(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv) +
				qQpipipd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv) +
				QpipiQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pxpxpypypd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv) +
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
				QpipiQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzppippipd(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, int num,
									  double D1deriv, double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv) +
				qQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv) +
				QpipiQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzpzpzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return qqpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) +
				qQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv) +
				QzzQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzsspd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								 double p1deriv, double p2deriv) {
		return -quzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzppippipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return -quzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				uzQpipipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzpzpzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, int num, double D1deriv,
								   double D2deriv, double p1deriv, double p2deriv) {
		return -quzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv, p2deriv) +
				uzQzzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ssspzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								 double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								 double p1deriv, double p2deriv) {
		return quzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppippispzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return quzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) -
				uzQpipipd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double pzpzspzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								   double p22, double D12, double D22, double R, int num, double D1deriv,
								   double D2deriv, double p1deriv, double p2deriv) {
		return quzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv) -
				uzQzzpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
						p2deriv);
	}

	public static double sppisppipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return upiupipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double spzspzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
								  double p22, double D12, double D22, double R, int num, double D1deriv, double D2deriv,
								  double p1deriv, double p2deriv) {
		return uzuzpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double sppippipzpd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return upiQpizpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double ppipzsppipd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriv,
									 double D2deriv, double p1deriv, double p2deriv) {
		return -upiQpizpd(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriv, D2deriv, p1deriv,
				p2deriv);
	}

	public static double ppipzppipzpd(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, int num, double D1deriv,
									  double D2deriv, double p1deriv, double p2deriv) {
		return QpizQpizpd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv, p2deriv);
	}

	public static double pxpypxpypd(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriv,
									double D2deriv, double p1deriv, double p2deriv) {
		return 0.5 *
				(ppippippippipd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
						p2deriv) -
						pxpxpypypd(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriv, D2deriv, p1deriv,
								p2deriv));
	}

	public static int f(int num) {
		if (num == 0 || num == 1) return 1 - num;
		else if (num == 2) return num;
		return -1;
	}

	public static double qqcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, double D11deriv, double D21deriv,
									 double p11deriv, double p21deriv, double D12deriv, double D22deriv,
									 double p12deriv, double p22deriv) {
		return 0;
	}

	public static double quzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, double D11deriv,
									  double D21deriv, double p11deriv, double p21deriv, double D12deriv,
									  double D22deriv, double p12deriv, double p22deriv) {
		return 0;
	}

	public static double qQpipicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, double D11deriv,
										 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										 double D22deriv, double p12deriv, double p22deriv) {
		return 0;
	}

	public static double qQzzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, double D11deriv,
									   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
									   double D22deriv, double p12deriv, double p22deriv) {
		return 0;
	}

	public static double h1(int a1, int a2, int b1, int b2, double Dn1deriv, double Dm2deriv, double pn1deriv,
							 double pm2deriv) {
		return (a1 * a2 + b1 * b2) * Dn1deriv * Dm2deriv + pn1deriv * pm2deriv;
	}

	public static double f1crossp2d(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2,
									 double R, double Dn1deriv, double Dm2deriv, double pn1deriv, double pm2deriv) {
		double f1 = f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);
		return 3 * Math.pow(f1, 5) * g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dn1deriv, pn1deriv, 0) *
				g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dm2deriv, pm2deriv, 1) -
				h1(a1, a2, b1, b2, Dn1deriv, Dm2deriv, pn1deriv, pm2deriv) * Math.pow(f1, 3);
	}

	public static double upiupicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, double D11deriv,
										 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										 double D22deriv, double p12deriv, double p22deriv) {
		return 0.5 * (+f1crossp2d(0, 0, 1, -1, D11, D12, p11, p12, R, D11deriv, D12deriv, p11deriv, p12deriv) -
				f1crossp2d(0, 0, 1, +1, D11, D12, p11, p12, R, D11deriv, D12deriv, p11deriv, p12deriv));
	}

	public static double uzuzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, double D11deriv,
									   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
									   double D22deriv, double p12deriv, double p22deriv) {
		return 0.25 * (+f1crossp2d(+1, -1, 0, 0, D11, D12, p11, p12, R, D11deriv, D12deriv, p11deriv, p12deriv) +
				f1crossp2d(-1, +1, 0, 0, D11, D12, p11, p12, R, D11deriv, D12deriv, p11deriv, p12deriv) -
				f1crossp2d(-1, -1, 0, 0, D11, D12, p11, p12, R, D11deriv, D12deriv, p11deriv, p12deriv) -
				f1crossp2d(+1, +1, 0, 0, D11, D12, p11, p12, R, D11deriv, D12deriv, p11deriv, p12deriv));
	}

	public static double upiQpizcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double R, double D11deriv,
										  double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										  double D22deriv, double p12deriv, double p22deriv) {
		return 0.25 * (+f1crossp2d(0, -1, +1, +1, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) -
				f1crossp2d(0, -1, +1, -1, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) +
				f1crossp2d(0, +1, +1, -1, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) -
				f1crossp2d(0, +1, +1, +1, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv));
	}

	public static double uzQpipicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double R, double D11deriv,
										  double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										  double D22deriv, double p12deriv, double p22deriv) {
		return 0.25 * (-f1crossp2d(+1, 0, 0, 2, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) +
				f1crossp2d(-1, 0, 0, 2, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) +
				f1crossp2d(+1, 0, 0, 0, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) -
				f1crossp2d(-1, 0, 0, 0, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv));
	}

	public static double uzQzzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, double D11deriv,
										double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										double D22deriv, double p12deriv, double p22deriv) {
		return 0.125 * (-f1crossp2d(+1, -2, 0, 0, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) +
				f1crossp2d(-1, -2, 0, 0, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) -
				f1crossp2d(+1, +2, 0, 0, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) +
				f1crossp2d(-1, +2, 0, 0, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) +
				2 * f1crossp2d(+1, 0, 0, 0, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv) -
				2 * f1crossp2d(-1, 0, 0, 0, D11, D22, p11, p22, R, D11deriv, D22deriv, p11deriv, p22deriv));
	}

	public static double QpipiQpipicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											 double p12, double p22, double D12, double D22, double R, double D11deriv,
											 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
											 double D22deriv, double p12deriv, double p22deriv) {
		return 0.125 * (+f1crossp2d(0, 0, +2, -2, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				f1crossp2d(0, 0, +2, +2, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				2 * f1crossp2d(0, 0, +2, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				2 * f1crossp2d(0, 0, 0, +2, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				2 * f1crossp2d(0, 0, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv));
	}

	public static double h2(int a1, int a2, double D21deriv, double D22deriv, double p21deriv, double p22deriv) {
		return a1 * a2 * D21deriv * D22deriv + p21deriv * p22deriv;
	}

	public static double f2crossp2d(int a1, int a2, int b1, int b2, double D21, double D22, double p21, double p22,
									 double R, double D21deriv, double D22deriv, double p21deriv, double p22deriv) {
		double f2 = f2(a1, a2, b1, b2, D21, D22, p21, p22, R);
		return 3 * Math.pow(f2, 5) * g2(a1, a2, b1, b2, D21, D22, p21, p22, R, D21deriv, p21deriv, 0) *
				g2(a1, a2, b1, b2, D21, D22, p21, p22, R, D22deriv, p22deriv, 1) -
				h2(a1, a2, D21deriv, D22deriv, p21deriv, p22deriv) * Math.pow(f2, 3);
	}

	public static double QxxQyycrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, double D11deriv,
										 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										 double D22deriv, double p12deriv, double p22deriv) {
		return 0.25 * (+f2crossp2d(0, 0, 2, 2, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				f2crossp2d(0, 0, 2, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				f2crossp2d(0, 0, 0, 2, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				f2crossp2d(0, 0, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv));
	}

	public static double QpipiQzzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, double D11deriv,
										   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										   double D22deriv, double p12deriv, double p22deriv) {
		return 0.125 * (+f1crossp2d(0, -2, +2, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				f1crossp2d(0, +2, +2, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				f1crossp2d(0, -2, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				f1crossp2d(0, +2, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				2 * f1crossp2d(0, 0, 2, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				2 * f1crossp2d(0, 0, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv));
	}

	public static double QzzQzzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, double D11deriv,
										 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										 double D22deriv, double p12deriv, double p22deriv) {
		return 0.0625 * (+f1crossp2d(+2, -2, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				f1crossp2d(+2, +2, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				f1crossp2d(-2, -2, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				f1crossp2d(-2, +2, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				2 * f1crossp2d(+2, 0, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				2 * f1crossp2d(-2, 0, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				2 * f1crossp2d(0, +2, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				2 * f1crossp2d(0, -2, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				4 * f1crossp2d(0, 0, 0, 0, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv));
	}

	public static double QpizQpizcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, double D11deriv,
										   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										   double D22deriv, double p12deriv, double p22deriv) {
		return 0.125 * (+f1crossp2d(+1, -1, +1, -1, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				f1crossp2d(+1, -1, +1, +1, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				f1crossp2d(+1, +1, +1, -1, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				f1crossp2d(+1, +1, +1, +1, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				f1crossp2d(-1, -1, +1, -1, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				f1crossp2d(-1, -1, +1, +1, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) +
				f1crossp2d(-1, +1, +1, -1, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv) -
				f1crossp2d(-1, +1, +1, +1, D21, D22, p21, p22, R, D21deriv, D22deriv, p21deriv, p22deriv));
	}

	public static double sssscrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, double D11deriv,
									   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
									   double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double ssppippicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, double D11deriv,
										   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										   double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) +
				qQpipicrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double sspzpzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, double D11deriv,
										 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										 double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) +
				qQzzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double ppippisscrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, double D11deriv,
										   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										   double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) +
				qQpipicrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv);
	}

	public static double pzpzsscrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, double D11deriv,
										 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										 double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) +
				qQzzcrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv);
	}

	public static double ppippippippicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											   double p12, double p22, double D12, double D22, double R,
											   double D11deriv, double D21deriv, double p11deriv, double p21deriv,
											   double D12deriv, double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) +
				qQpipicrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv) +
				qQpipicrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv) +
				QpipiQpipicrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double pxpxpypycrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, double D11deriv,
										   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										   double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) +
				qQpipicrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv) +
				qQpipicrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv) +
				QxxQyycrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double ppippipzpzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											 double p12, double p22, double D12, double D22, double R, double D11deriv,
											 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
											 double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) +
				qQzzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv) +
				qQpipicrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv) +
				QpipiQzzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double pzpzppippicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											 double p12, double p22, double D12, double D22, double R, double D11deriv,
											 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
											 double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) +
				qQpipicrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv) +
				qQzzcrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv) +
				QpipiQzzcrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv);
	}

	public static double pzpzpzpzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, double D11deriv,
										   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										   double D22deriv, double p12deriv, double p22deriv) {
		return qqcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) +
				qQzzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv) +
				qQzzcrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv) +
				QzzQzzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double spzsscrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, double D11deriv,
										double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										double D22deriv, double p12deriv, double p22deriv) {
		return -quzcrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
				p22deriv,
				D11deriv, D21deriv, p11deriv, p21deriv);
	}

	public static double spzppippicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, double D11deriv,
											double D21deriv, double p11deriv, double p21deriv, double D12deriv,
											double D22deriv, double p12deriv, double p22deriv) {
		return -quzcrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
				p22deriv,
				D11deriv, D21deriv, p11deriv, p21deriv) +
				uzQpipicrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double spzpzpzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double R, double D11deriv,
										  double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										  double D22deriv, double p12deriv, double p22deriv) {
		return -quzcrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
				p22deriv,
				D11deriv, D21deriv, p11deriv, p21deriv) +
				uzQzzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double ssspzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, double D11deriv,
										double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										double D22deriv, double p12deriv, double p22deriv) {
		return quzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double ppippispzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, double D11deriv,
											double D21deriv, double p11deriv, double p21deriv, double D12deriv,
											double D22deriv, double p12deriv, double p22deriv) {
		return quzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) -
				uzQpipicrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv);
	}

	public static double pzpzspzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double R, double D11deriv,
										  double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										  double D22deriv, double p12deriv, double p22deriv) {
		return quzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv, p21deriv,
				D12deriv, D22deriv, p12deriv, p22deriv) -
				uzQzzcrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
						p22deriv, D11deriv, D21deriv, p11deriv, p21deriv);
	}

	public static double sppisppicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, double D11deriv,
										   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										   double D22deriv, double p12deriv, double p22deriv) {
		return upiupicrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
				p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double spzspzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, double D11deriv,
										 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										 double D22deriv, double p12deriv, double p22deriv) {
		return uzuzcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
				p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double sppippipzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, double D11deriv,
											double D21deriv, double p11deriv, double p21deriv, double D12deriv,
											double D22deriv, double p12deriv, double p22deriv) {
		return upiQpizcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
				p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double ppipzsppicrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, double D11deriv,
											double D21deriv, double p11deriv, double p21deriv, double D12deriv,
											double D22deriv, double p12deriv, double p22deriv) {
		return -upiQpizcrossp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, D12deriv, D22deriv, p12deriv,
				p22deriv, D11deriv, D21deriv, p11deriv, p21deriv);
	}

	public static double ppipzppipzcrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											 double p12, double p22, double D12, double D22, double R, double D11deriv,
											 double D21deriv, double p11deriv, double p21deriv, double D12deriv,
											 double D22deriv, double p12deriv, double p22deriv) {
		return QpizQpizcrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
				p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
	}

	public static double pxpypxpycrossp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, double D11deriv,
										   double D21deriv, double p11deriv, double p21deriv, double D12deriv,
										   double D22deriv, double p12deriv, double p22deriv) {
		return 0.5 *
				(ppippippippicrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv, p11deriv,
						p21deriv, D12deriv, D22deriv, p12deriv, p22deriv) -
						pxpxpypycrossp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, D11deriv, D21deriv,
								p11deriv, p21deriv, D12deriv, D22deriv, p12deriv, p22deriv));
	}

	public static double h0(int a2, int b2, double Dn2, double p01, double pn2, double R, double Dn2deriva,
							 double pn2deriva, double Dn2derivb, double pn2derivb, double Dn2deriv2,
							 double pn2deriv2) {
		return a2 * Dn2deriv2 * (R + a2 * Dn2) + a2 * a2 * Dn2deriva * Dn2derivb +
				b2 * b2 * (Dn2 * Dn2deriv2 + Dn2deriva * Dn2derivb) + (p01 + pn2) * pn2deriv2 + pn2deriva * pn2derivb;
	}

	public static double f0diagp2dpartial(int a2, int b2, double Dn2, double p01, double pn2, double R,
										   double Dn2deriva, double pn2deriva, double Dn2derivb, double pn2derivb,
										   double Dn2deriv2, double pn2deriv2) {
		double f0 = f0(a2, b2, Dn2, p01, pn2, R);

		return 3 * Math.pow(f0, 5) * g0(a2, b2, Dn2, p01, pn2, R, Dn2deriva, pn2deriva) *
				g0(a2, b2, Dn2, p01, pn2, R, Dn2derivb, pn2derivb)
				- h0(a2, b2, Dn2, p01, pn2, R, Dn2deriva, pn2deriva, Dn2derivb, pn2derivb, Dn2deriv2, pn2deriv2) *
				Math.pow(f0, 3);
	}

	public static double f0diagp2d(int a2, int b2, double Dn2, double p01, double pn2, double R, double Dn2deriva,
									double pn2deriva, double Dn2derivb, double pn2derivb, double Dn2deriv2,
									double pn2deriv2, int type) {
		switch (type) {
			case 0:
				return 0;
			case 1:
			case 2:
				return f0diagp2dpartial(a2, b2, Dn2, p01, pn2, R, Dn2deriva, pn2deriva, Dn2derivb, pn2derivb,
						Dn2deriv2,
						pn2deriv2);
		}
		return 0;
	}


	public static double qqdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02, double p12,
									double p22, double D12, double D22, double R, int num, double D1deriva,
									double D2deriva, double p1deriva, double p2deriva, double D1derivb,
									double D2derivb,
									double p1derivb, double p2derivb, double D1deriv2, double D2deriv2,
									double p1deriv2,
									double p2deriv2) {
		return 0;
	}

	public static double quzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
									 double p12,
									 double p22, double D12, double D22, double R, int num, double D1deriva,
									 double D2deriva, double p1deriva, double p2deriva, double D1derivb,
									 double D2derivb, double p1derivb, double p2derivb, double D1deriv2,
									 double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.5 * (
				+f0diagp2d(+1, 0, D12, p01, p12, R, D1deriva, p1deriva, D1derivb, p1derivb, D1deriv2, p1deriv2, num)
						- f0diagp2d(-1, 0, D12, p01, p12, R, D1deriva, p1deriva, D1derivb, p1derivb, D1deriv2,
						p1deriv2,
						num)
		);
	}

	public static double qQpipidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.5 * (
				+f0diagp2d(0, 2, D22, p01, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						- f0diagp2d(0, 0, D22, p01, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2,
						num)
		);
	}

	public static double qQzzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, int num,
									  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
									  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
									  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.25 * (
				+f0diagp2d(+2, 0, D22, p01, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						+ f0diagp2d(-2, 0, D22, p01, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
						p2deriv2,
						num)
						- 2 *
						f0diagp2d(0, 0, D22, p01, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2,
								num)

		);
	}

	public static double h1(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R,
							 double Dxderiva, double pxderiva, double Dxderivb, double pxderivb, double Dxderiv2,
							 double pxderiv2, int type) {
		switch (type) {
			case 0:
				return a1 * a1 * Dxderiva * Dxderivb + a1 * Dxderiv2 * (R + a1 * Dn1 + a2 * Dm2) +
						b1 * b1 * Dxderiva * Dxderivb + b1 * Dxderiv2 * (b1 * Dn1 + b2 * Dm2) + pxderiva * pxderivb +
						(pn1 + pm2) * pxderiv2;
			case 1:
				return a2 * a2 * Dxderiva * Dxderivb + a2 * Dxderiv2 * (R + a1 * Dn1 + a2 * Dm2) +
						b2 * b2 * Dxderiva * Dxderivb + b2 * Dxderiv2 * (b1 * Dn1 + b2 * Dm2) + pxderiva * pxderivb +
						(pn1 + pm2) * pxderiv2;

		}

		return 0;
	}

	public static double f1diagp2dpartial(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1,
										   double pm2, double R, double Dxderiva, double pxderiva, double Dxderivb,
										   double pxderivb, double Dxderiv2, double pxderiv2, int type) {
		double f1 = f1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);

		return 3 * Math.pow(f1, 5) * g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, type) *
				g1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderivb, pxderivb, type)
				- Math.pow(f1, 3) *
				h1(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, Dxderivb, pxderivb, Dxderiv2, pxderiv2,
						type);
	}

	public static double f1diagp2d(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2,
									double R, double Dxderiva, double pxderiva, double Dxderivb, double pxderivb,
									double Dxderiv2, double pxderiv2, int type) {
		switch (type) {
			case 0:
			case 1:
				return f1diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, Dxderivb, pxderivb,
						Dxderiv2, pxderiv2, type);
			case 2:
				return f1diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, Dxderivb, pxderivb,
						Dxderiv2, pxderiv2, 0)
						+
						f1diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, Dxderivb, pxderivb,
								Dxderiv2, pxderiv2, 1)
						+ f1crossp2d(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, Dxderivb, pxderiva, pxderivb) * 2;
		}

		return 0;
	}

	public static double upiupidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.5 * (
				+f1diagp2d(0, 0, 1, -1, D11, D12, p11, p12, R, D1deriva, p1deriva, D1derivb, p1derivb, D1deriv2,
						p1deriv2, num)
						-
						f1diagp2d(0, 0, 1, +1, D11, D12, p11, p12, R, D1deriva, p1deriva, D1derivb, p1derivb, D1deriv2,
								p1deriv2, num)
		);
	}

	public static double uzuzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
									  double p12, double p22, double D12, double D22, double R, int num,
									  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
									  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
									  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.25 * (
				+f1diagp2d(+1, -1, 0, 0, D11, D12, p11, p12, R, D1deriva, p1deriva, D1derivb, p1derivb, D1deriv2,
						p1deriv2, num)
						+
						f1diagp2d(-1, +1, 0, 0, D11, D12, p11, p12, R, D1deriva, p1deriva, D1derivb, p1derivb,
								D1deriv2,
								p1deriv2, num)
						-
						f1diagp2d(-1, -1, 0, 0, D11, D12, p11, p12, R, D1deriva, p1deriva, D1derivb, p1derivb,
								D1deriv2,
								p1deriv2, num)
						-
						f1diagp2d(+1, +1, 0, 0, D11, D12, p11, p12, R, D1deriva, p1deriva, D1derivb, p1derivb,
								D1deriv2,
								p1deriv2, num)
		);
	}

	public static double varf1diagp2d(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2,
									   double R, double Dn1deriva, double pn1deriva, double Dn1derivb,
									   double pn1derivb,
									   double Dn1deriv2, double pn1deriv2, double Dm2deriva, double pm2deriva,
									   double Dm2derivb, double pm2derivb, double Dm2deriv2, double pm2deriv2,
									   int type) {
		switch (type) {
			case 0:
				return f1diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dn1deriva, pn1deriva, Dn1derivb,
						pn1derivb, Dn1deriv2, pn1deriv2, type);
			case 1:
				return f1diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dm2deriva, pm2deriva, Dm2derivb,
						pm2derivb, Dm2deriv2, pm2deriv2, type);
			case 2:
				return f1diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dn1deriva, pn1deriva, Dn1derivb,
						pn1derivb, Dn1deriv2, pn1deriv2, 0)
						+ f1diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dm2deriva, pm2deriva, Dm2derivb,
						pm2derivb, Dm2deriv2, pm2deriv2, 1)
						+ f1crossp2d(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dn1deriva, Dm2derivb, pn1deriva, pm2derivb)
						+ f1crossp2d(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dn1derivb, Dm2deriva, pn1derivb,
						pm2deriva);
		}

		return 0;
	}

	public static double upiQpizdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										 double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										 double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.25 * (
				+varf1diagp2d(0, -1, +1, +1, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb, D1deriv2,
						p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						- varf1diagp2d(0, -1, +1, -1, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						+ varf1diagp2d(0, +1, +1, -1, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						- varf1diagp2d(0, +1, +1, +1, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
		);
	}

	public static double uzQpipidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										 double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										 double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.25 * (
				-varf1diagp2d(+1, 0, 0, 2, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb, D1deriv2,
						p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						+ varf1diagp2d(-1, 0, 0, 2, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						+ varf1diagp2d(+1, 0, 0, 0, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						- varf1diagp2d(-1, 0, 0, 0, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
		);
	}

	public static double uzQzzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
									   double p12, double p22, double D12, double D22, double R, int num,
									   double D1deriva, double D2deriva, double p1deriva, double p2deriva,
									   double D1derivb, double D2derivb, double p1derivb, double p2derivb,
									   double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.125 * (
				-varf1diagp2d(+1, -2, 0, 0, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb, D1deriv2,
						p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						+ varf1diagp2d(-1, -2, 0, 0, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						- varf1diagp2d(+1, +2, 0, 0, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						+ varf1diagp2d(-1, +2, 0, 0, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						+ 2 * varf1diagp2d(+1, 0, 0, 0, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
						- 2 * varf1diagp2d(-1, 0, 0, 0, D11, D22, p11, p22, R, D1deriva, p1deriva, D1derivb, p1derivb,
						D1deriv2, p1deriv2, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2, p2deriv2, num)
		);
	}


	public static double QpipiQpipidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, int num,
											double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.125 * (
				+f1diagp2d(0, 0, +2, -2, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
						p2deriv2, num)
						+
						f1diagp2d(0, 0, +2, +2, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
								D2deriv2,
								p2deriv2, num)
						- 2 *
						f1diagp2d(0, 0, +2, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
						- 2 *
						f1diagp2d(0, 0, 0, +2, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
						+ 2 *
						f1diagp2d(0, 0, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
		);
	}

	public static double h2(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2, double R,
							 double Dxderiva, double pxderiva, double Dxderivb, double pxderivb, double Dxderiv2,
							 double pxderiv2, int type) {
		switch (type) {
			case 0:
				return a1 * a1 * Dxderiva * Dxderivb + a1 * Dxderiv2 * (R + a1 * Dn1 + a2 * Dm2) +
						b1 * b1 * Dxderiva * Dxderivb + b1 * Dxderiv2 * b1 * Dn1 + pxderiva * pxderivb +
						(pn1 + pm2) * pxderiv2;
			case 1:
				return a2 * a2 * Dxderiva * Dxderivb + a2 * Dxderiv2 * (R + a1 * Dn1 + a2 * Dm2) +
						b2 * b2 * Dxderiva * Dxderivb + b2 * Dxderiv2 * b2 * Dm2 + pxderiva * pxderivb +
						(pn1 + pm2) * pxderiv2;

		}

		return 0;
	}

	public static double f2diagp2dpartial(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1,
										   double pm2, double R, double Dxderiva, double pxderiva, double Dxderivb,
										   double pxderivb, double Dxderiv2, double pxderiv2, int type) {
		double f2 = f2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R);

		return 3 * Math.pow(f2, 5) * g2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, type) *
				g2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderivb, pxderivb, type)
				- Math.pow(f2, 3) *
				h2(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, Dxderivb, pxderivb, Dxderiv2, pxderiv2,
						type);
	}

	public static double f2diagp2d(int a1, int a2, int b1, int b2, double Dn1, double Dm2, double pn1, double pm2,
									double R, double Dxderiva, double pxderiva, double Dxderivb, double pxderivb,
									double Dxderiv2, double pxderiv2, int type) {
		switch (type) {
			case 0:
			case 1:
				return f2diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, Dxderivb, pxderivb,
						Dxderiv2, pxderiv2, type);
			case 2:
				return f2diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, Dxderivb, pxderivb,
						Dxderiv2, pxderiv2, 0)
						+
						f2diagp2dpartial(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, pxderiva, Dxderivb, pxderivb,
								Dxderiv2, pxderiv2, 1)
						+ f2crossp2d(a1, a2, b1, b2, Dn1, Dm2, pn1, pm2, R, Dxderiva, Dxderivb, pxderiva, pxderivb) * 2;
		}

		return 0;
	}


	public static double QxxQyydiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.25 * (
				+f2diagp2d(0, 0, 2, 2, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
						p2deriv2, num)
						- f2diagp2d(0, 0, 2, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2,
						p2deriv2, num)
						- f2diagp2d(0, 0, 0, 2, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2,
						p2deriv2, num)
						+ f2diagp2d(0, 0, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2,
						p2deriv2, num)
		);
	}

	public static double QpipiQzzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double R, int num,
										  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.125 * (
				+f1diagp2d(0, -2, +2, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
						p2deriv2, num)
						+
						f1diagp2d(0, +2, +2, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
								D2deriv2,
								p2deriv2, num)
						-
						f1diagp2d(0, -2, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
						-
						f1diagp2d(0, +2, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
						- 2 *
						f1diagp2d(0, 0, 2, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
						+ 2 *
						f1diagp2d(0, 0, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
		);

	}

	public static double QzzQzzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.0625 * (
				+f1diagp2d(+2, -2, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
						p2deriv2, num)
						+
						f1diagp2d(+2, +2, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
								D2deriv2,
								p2deriv2, num)
						+
						f1diagp2d(-2, -2, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
								D2deriv2,
								p2deriv2, num)
						+
						f1diagp2d(-2, +2, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
								D2deriv2,
								p2deriv2, num)
						- 2 *
						f1diagp2d(+2, 0, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
						- 2 *
						f1diagp2d(-2, 0, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
						- 2 *
						f1diagp2d(0, +2, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
						- 2 *
						f1diagp2d(0, -2, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
						+ 4 *
						f1diagp2d(0, 0, 0, 0, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
								p2deriv2, num)
		);
	}

	public static double QpizQpizdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double R, int num,
										  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.125 * (
				+f1diagp2d(+1, -1, +1, -1, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb, D2deriv2,
						p2deriv2, num)
						- f1diagp2d(+1, -1, +1, +1, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2, p2deriv2, num)
						- f1diagp2d(+1, +1, +1, -1, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2, p2deriv2, num)
						+ f1diagp2d(+1, +1, +1, +1, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2, p2deriv2, num)
						- f1diagp2d(-1, -1, +1, -1, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2, p2deriv2, num)
						+ f1diagp2d(-1, -1, +1, +1, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2, p2deriv2, num)
						+ f1diagp2d(-1, +1, +1, -1, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2, p2deriv2, num)
						- f1diagp2d(-1, +1, +1, +1, D21, D22, p21, p22, R, D2deriva, p2deriva, D2derivb, p2derivb,
						D2deriv2, p2deriv2, num)
		);
	}

	public static double ssssdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										double p12, double p22, double D12, double D22, double R, int num,
										double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double ssppippidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, int num,
											double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQpipidiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double sspzpzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double R, int num,
										  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQzzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double ppippissdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, int num,
											double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQpipidiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva,
						p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double pzpzssdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double R, int num,
										  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQzzdiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double ppippippippidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
												double p12, double p22, double D12, double D22, double R, int num,
												double D1deriva, double D2deriva, double p1deriva, double p2deriva,
												double D1derivb, double D2derivb, double p1derivb, double p2derivb,
												double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQpipidiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQpipidiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva,
						p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				QpipiQpipidiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva,
						p1deriva, p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2,
						p2deriv2);
	}

	public static double pxpxpypydiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, int num,
											double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQpipidiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQpipidiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva,
						p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				QxxQyydiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double ppippipzpzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											  double p12, double p22, double D12, double D22, double R, int num,
											  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQzzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQpipidiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva,
						p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				QpipiQzzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double pzpzppippidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											  double p12, double p22, double D12, double D22, double R, int num,
											  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQpipidiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQzzdiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				QpipiQzzdiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva,
						p1deriva, p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2,
						p2deriv2);
	}

	public static double pzpzpzpzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, int num,
											double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return qqdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQzzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				qQzzdiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				QzzQzzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double spzssdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										 double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										 double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return -quzdiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double spzppippidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											 double p12, double p22, double D12, double D22, double R, int num,
											 double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											 double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											 double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return -quzdiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				uzQpipidiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double spzpzpzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, int num,
										   double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										   double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										   double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return -quzdiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) +
				uzQzzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double ssspzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										 double p12, double p22, double D12, double D22, double R, int num,
										 double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										 double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										 double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return quzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double ppippispzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											 double p12, double p22, double D12, double D22, double R, int num,
											 double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											 double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											 double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return quzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) -
				uzQpipidiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva,
						p1deriva, p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2,
						p2deriv2);
	}

	public static double pzpzspzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										   double p12, double p22, double D12, double D22, double R, int num,
										   double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										   double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										   double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return quzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2) -
				uzQzzdiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva, p1deriva,
						p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double sppisppidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, int num,
											double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return upiupidiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double spzspzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
										  double p12, double p22, double D12, double D22, double R, int num,
										  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
										  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
										  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return uzuzdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double sppippipzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											 double p12, double p22, double D12, double D22, double R, int num,
											 double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											 double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											 double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return upiQpizdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double ppipzsppidiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											 double p12, double p22, double D12, double D22, double R, int num,
											 double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											 double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											 double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return -upiQpizdiagp2d(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, R, f(num), D1deriva, D2deriva,
				p1deriva, p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double ppipzppipzdiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											  double p12, double p22, double D12, double D22, double R, int num,
											  double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											  double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											  double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return QpizQpizdiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva, p1deriva,
				p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2, p2deriv2);
	}

	public static double pxpypxpydiagp2d(double p01, double p11, double p21, double D11, double D21, double p02,
											double p12, double p22, double D12, double D22, double R, int num,
											double D1deriva, double D2deriva, double p1deriva, double p2deriva,
											double D1derivb, double D2derivb, double p1derivb, double p2derivb,
											double D1deriv2, double D2deriv2, double p1deriv2, double p2deriv2) {
		return 0.5 *
				(ppippippippidiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva,
						p1deriva, p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2, p1deriv2,
						p2deriv2) -
						pxpxpypydiagp2d(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, R, num, D1deriva, D2deriva,
								p1deriva, p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2,
								p1deriv2, p2deriv2));
	}
}
