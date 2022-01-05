package nddo.scf;

import nddo.defaults.NDDO6G;
import nddo.structs.OrbitalProperties;

public class STO6G extends LCGTO {
	private final static double[] exp1 =
			new double[]{6.51095E-2, 1.58088E-1, 4.07099E-1, 1.18506, 4.23592, 2.31030E+1};

	private final static double[] coeff1 =
			new double[]{1.30334E-1, 4.16492E-1, 3.70563E-1, 1.68538E-1, 4.93615E-2, 9.16360E-3};

	private final static double[] exp2 =
			new double[]{4.85690E-2, 1.05960E-1, 2.43977E-1, 6.34142E-1, 2.04036, 1.03087E+1};

	private final static double[] coeff2s =
			new double[]{2.40706E-1, 5.95117E-1, 2.50242E-1, -3.37854E-2, -4.69917E-2, -1.32528E-2};

	private final static double[] coeff2p =
			new double[]{1.01708E-1, 4.25860E-1, 4.18036E-1, 1.73897E-1, 3.76794E-2, 3.75970E-3};

	public final double zeta;

	public STO6G(OrbitalProperties op, double[] coordinates, double zeta) {
		super(op, coordinates, getExps(op.shell, zeta), getCoeffs(op.shell, op.L));
		this.zeta = zeta;
	}

	public STO6G(STO6G sto6G) {
		super(sto6G);
		this.zeta = sto6G.zeta;
	}

	private static double[] getExps(int shell, double zeta) {
		switch (shell) {
			case 1:
				double[] exp1 = new double[6];
				for (int i = 0; i < 6; i++) {
					exp1[i] = STO6G.exp1[i] * zeta * zeta;
				}
				return exp1;
			case 2:
				double[] exp2 = new double[6];
				for (int i = 0; i < 6; i++) {
					exp2[i] = STO6G.exp2[i] * zeta * zeta;
				}
				return exp2;
		}

		throw new IllegalArgumentException("Illegal atom");
	}

	private static double[] getCoeffs(int shell, int L) {
		switch (shell) {
			case 1:
				return coeff1.clone();
			case 2:
				switch (L) {
					case 0:
						return coeff2s.clone();
					case 1:
						return coeff2p.clone();
				}
		}

		throw new IllegalArgumentException("Illegal atom" + shell + L);
	}

	public static double Spd(STO6G X1, STO6G X2, int num, int type) {
		if (num == -1) {
			return 0;
		}

		int hasA = X1.getL() == type && num != 1 ? 1 : 0;
		int hasB = X2.getL() == type && num != 0 ? 1 : 0;

		hasB <<= 1;
		int alltogether = hasA + hasB - 1;

		double Sderiv = 0;

		for (int i = 0; i < X1.getn(); i++) {
			for (int j = 0; j < X2.getn(); j++) {

				switch (alltogether) {
					case 0:
						Sderiv += X1.gaussExponents[i] * 2 / X1.zeta * X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphapd(X1.getGaussArray()[i], X2.getGaussArray()[j], 0);
						break;
					case 1:
						Sderiv += X2.gaussExponents[j] * 2 / X2.zeta * X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphapd(X1.getGaussArray()[i], X2.getGaussArray()[j], 1);
						break;
					case 2:
						Sderiv += X1.gaussExponents[i] * 2 / X1.zeta * X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphapd(X1.getGaussArray()[i], X2.getGaussArray()[j], 0);

						Sderiv += X2.gaussExponents[j] * 2 / X2.zeta * X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphapd(X1.getGaussArray()[i], X2.getGaussArray()[j], 1);
				}
			}

		}

		return Sderiv * X1.getN() * X2.getN();
	}

	private static double Szetadiagp2d(STO6G X1, STO6G X2, int type) {
		double Sderiv = 0;

		for (int i = 0; i < X1.getn(); i++) {
			for (int j = 0; j < X2.getn(); j++) {

				switch (type) {
					case 0:
						Sderiv += X1.gaussExponents[i] * X1.gaussExponents[i] * 4 / (X1.zeta * X1.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphadiagp2d(X1.getGaussArray()[i], X2.getGaussArray()[j], 0);
						Sderiv += 2 * X1.gaussExponents[i] / (X1.zeta * X1.zeta) * X1.getCoeffArray()[i] *
								X2.getCoeffArray()[j] *
								GTO.Salphapd(X1.getGaussArray()[i], X2.getGaussArray()[j], 0);

						break;
					case 1:
						Sderiv += X2.gaussExponents[j] * X2.gaussExponents[j] * 4 / (X2.zeta * X2.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphadiagp2d(X1.getGaussArray()[i], X2.getGaussArray()[j], 1);
						Sderiv += 2 * X2.gaussExponents[j] / (X2.zeta * X2.zeta) * X1.getCoeffArray()[i] *
								X2.getCoeffArray()[j] *
								GTO.Salphapd(X1.getGaussArray()[i], X2.getGaussArray()[j], 1);

						break;
					case 2:
						Sderiv += X1.gaussExponents[i] * X1.gaussExponents[i] * 4 / (X1.zeta * X1.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphadiagp2d(X1.getGaussArray()[i], X2.getGaussArray()[j], 0);
						Sderiv += 2 * X1.gaussExponents[i] / (X1.zeta * X1.zeta) * X1.getCoeffArray()[i] *
								X2.getCoeffArray()[j] *
								GTO.Salphapd(X1.getGaussArray()[i], X2.getGaussArray()[j], 0);
						Sderiv += X2.gaussExponents[j] * X2.gaussExponents[j] * 4 / (X2.zeta * X2.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphadiagp2d(X1.getGaussArray()[i], X2.getGaussArray()[j], 1);
						Sderiv += 2 * X2.gaussExponents[j] / (X2.zeta * X2.zeta) * X1.getCoeffArray()[i] *
								X2.getCoeffArray()[j] *
								GTO.Salphapd(X1.getGaussArray()[i], X2.getGaussArray()[j], 1);
						Sderiv += 2 * X1.gaussExponents[i] * X2.gaussExponents[j] * 4 / (X1.zeta * X2.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphacrossp2d(X1.getGaussArray()[i], X2.getGaussArray()[j]);
				}
			}
		}

		return Sderiv * X1.getN() * X2.getN();
	}

	private static double Szetacrossp2d(STO6G X1, STO6G X2) {
		double Sderiv = 0;

		for (int i = 0; i < X1.getn(); i++) {
			for (int j = 0; j < X2.getn(); j++) {
				Sderiv += X1.gaussExponents[i] * X2.gaussExponents[j] * 4 / (X1.zeta * X2.zeta) *
						X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
						GTO.Salphacrossp2d(X1.getGaussArray()[i], X2.getGaussArray()[j]);
			}
		}

		return Sderiv * X1.getN() * X2.getN();
	}

	public static double Sp2d(STO6G a, STO6G b, int num1, int type1, int num2, int type2) {
		double returnval = 0;

		if (num1 != -1 && num2 != -1) {
			int hasA1 = a.getL() == type1 && num1 != 1 ? 1 : 0;
			int hasB1 = b.getL() == type1 && num1 != 0 ? 1 : 0;

			int hasA2 = a.getL() == type2 && num2 != 1 ? 1 : 0;
			int hasB2 = b.getL() == type2 && num2 != 0 ? 1 : 0;


			if (num1 == num2 && type1 == type2 && hasA1 + hasB1 > 0) {
				int alltogether = hasA1 + (hasB1 << 1) - 1;

				returnval = Szetadiagp2d(a, b, alltogether);
			}

			else if (hasA1 + hasB1 > 0 && hasA2 + hasB2 > 0) {

				if (((hasA1 ^ hasB1) & (hasA2 ^ hasB2)) == 1) {
					returnval = Szetacrossp2d(a, b);
				}
			}
		}

		return returnval;
	}

	public static double Spgd(STO6G a, STO6G b, int num, int type, int tau) {
		if (num == -1) {
			return 0;
		}

		int hasA = a.getL() == type && num != 1 ? 1 : 0;
		int hasB = b.getL() == type && num != 0 ? 1 : 0;

		hasB <<= 1;
		int alltogether = hasA + hasB - 1;

		double Sderiv = 0;

		for (int i = 0; i < a.getn(); i++) {
			for (int j = 0; j < b.getn(); j++) {
				switch (alltogether) {
					case 0:
						Sderiv += a.gaussExponents[i] * 2 / a.zeta * a.getCoeffArray()[i] * b.getCoeffArray()[j] *
								GTO.Salphapgd(a.getGaussArray()[i], b.getGaussArray()[j], 0, tau);
						break;
					case 1:
						Sderiv += b.gaussExponents[j] * 2 / b.zeta * a.getCoeffArray()[i] * b.getCoeffArray()[j] *
								GTO.Salphapgd(a.getGaussArray()[i], b.getGaussArray()[j], 1, tau);
						break;
					case 2:
						Sderiv += a.gaussExponents[i] * 2 / a.zeta * a.getCoeffArray()[i] * b.getCoeffArray()[j] *
								GTO.Salphapgd(a.getGaussArray()[i], b.getGaussArray()[j], 0, tau);

						Sderiv += b.gaussExponents[j] * 2 / b.zeta * a.getCoeffArray()[i] * b.getCoeffArray()[j] *
								GTO.Salphapgd(a.getGaussArray()[i], b.getGaussArray()[j], 1, tau);
				}
			}

		}

		return Sderiv * a.getN() * b.getN();
	}

	private static double Scrossp2gd(STO6G X1, STO6G X2, int tau) {
		double Sderiv = 0;

		for (int i = 0; i < X1.getn(); i++) {
			for (int j = 0; j < X2.getn(); j++) {
				Sderiv +=
						X1.gaussExponents[i] * X2.gaussExponents[j] * 4 / (X1.zeta * X2.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphacrossp2gd(X1.getGaussArray()[i], X2.getGaussArray()[j], tau);
			}
		}

		return Sderiv * X1.getN() * X2.getN();
	}

	private static double Sdiagp2gd(STO6G X1, STO6G X2, int type, int tau) {
		double Sderiv = 0;

		for (int i = 0; i < X1.getn(); i++) {
			for (int j = 0; j < X2.getn(); j++) {

				switch (type) {
					case 0:
						Sderiv += X1.gaussExponents[i] * X1.gaussExponents[i] * 4 / (X1.zeta * X1.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphadiagp2gd(X1.getGaussArray()[i], X2.getGaussArray()[j], 0, tau);
						Sderiv += 2 * X1.gaussExponents[i] / (X1.zeta * X1.zeta) * X1.getCoeffArray()[i] *
								X2.getCoeffArray()[j] *
								GTO.Salphapgd(X1.getGaussArray()[i], X2.getGaussArray()[j], 0, tau);

						break;
					case 1:
						Sderiv += X2.gaussExponents[j] * X2.gaussExponents[j] * 4 / (X2.zeta * X2.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphadiagp2gd(X1.getGaussArray()[i], X2.getGaussArray()[j], 1, tau);
						Sderiv += 2 * X2.gaussExponents[j] / (X2.zeta * X2.zeta) * X1.getCoeffArray()[i] *
								X2.getCoeffArray()[j] *
								GTO.Salphapgd(X1.getGaussArray()[i], X2.getGaussArray()[j], 1, tau);

						break;
					case 2:
						Sderiv += X1.gaussExponents[i] * X1.gaussExponents[i] * 4 / (X1.zeta * X1.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphadiagp2gd(X1.getGaussArray()[i], X2.getGaussArray()[j], 0, tau);
						Sderiv += 2 * X1.gaussExponents[i] / (X1.zeta * X1.zeta) * X1.getCoeffArray()[i] *
								X2.getCoeffArray()[j] *
								GTO.Salphapgd(X1.getGaussArray()[i], X2.getGaussArray()[j], 0, tau);
						Sderiv += X2.gaussExponents[j] * X2.gaussExponents[j] * 4 / (X2.zeta * X2.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphadiagp2gd(X1.getGaussArray()[i], X2.getGaussArray()[j], 1, tau);
						Sderiv += 2 * X2.gaussExponents[j] / (X2.zeta * X2.zeta) * X1.getCoeffArray()[i] *
								X2.getCoeffArray()[j] *
								GTO.Salphapgd(X1.getGaussArray()[i], X2.getGaussArray()[j], 1, tau);
						Sderiv += 2 * X1.gaussExponents[i] * X2.gaussExponents[j] * 4 / (X1.zeta * X2.zeta) *
								X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
								GTO.Salphacrossp2gd(X1.getGaussArray()[i], X2.getGaussArray()[j], tau);


				}
			}

		}

		return Sderiv * X1.getN() * X2.getN();
	}

	public static double Sp2gd(NDDO6G a, NDDO6G b, int num1, int type1, int num2, int type2, int tau) {
		double returnval = 0;

		if (num1 != -1 && num2 != -1) {
			int hasA1 = a.getL() == type1 && num1 != 1 ? 1 : 0;
			int hasB1 = b.getL() == type1 && num1 != 0 ? 1 : 0;

			int hasA2 = a.getL() == type2 && num2 != 1 ? 1 : 0;
			int hasB2 = b.getL() == type2 && num2 != 0 ? 1 : 0;


			if (num1 == num2 && type1 == type2 && hasA1 + hasB1 > 0) {
				int alltogether = hasA1 + (hasB1 << 1) - 1;

				returnval = STO6G.Sdiagp2gd(a, b, alltogether, tau);
			}

			else if (hasA1 + hasB1 > 0 && hasA2 + hasB2 > 0) {

				if (((hasA1 ^ hasB1) & (hasA2 ^ hasB2)) == 1) {
					returnval = STO6G.Scrossp2gd(a, b, tau);
				}
			}
		}

		return returnval;
	}
}
