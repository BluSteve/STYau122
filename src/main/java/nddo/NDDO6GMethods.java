package nddo;

import nddo.math.ERI;
import nddo.scf.LCGTO;
import nddo.scf.STO6G;

public class NDDO6GMethods implements NDDOOrbitalMethods<NDDO6G> {
	@Override
	public double OneCenterERI(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d) {
		return ERI.OneCenterERI(a, b, c, d);
	}

	@Override
	public double G(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d) {
		double[] coeffA = a.decomposition(a.coordinates, c.coordinates);
		double[] coeffB = b.decomposition(a.coordinates, c.coordinates);
		double[] coeffC = c.decomposition(a.coordinates, c.coordinates);
		double[] coeffD = d.decomposition(a.coordinates, c.coordinates);

		double sum2 = 0;

		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		for (int i = 0; i < coeffA.length; i++)
			for (int j = 0; j < coeffB.length; j++)
				for (int k = 0; k < coeffC.length; k++)
					for (int l = 0; l < coeffD.length; l++)
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum2 += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
									ERI.LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}

		return sum2 * Constants.eV;
	}

	@Override
	public double Ggd(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau) {
		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double[] coeffAderiv = a.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffBderiv = b.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffCderiv = c.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffDderiv = d.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);

		double sum = 0;

		if (Math.abs(a.getCoords()[0] - c.getCoords()[0]) < 1E-3 &&
				Math.abs(a.getCoords()[1] - c.getCoords()[1]) < 1E-3) {
			coeffA = a.decompositionvar(a.getCoords(), c.getCoords());
			coeffB = b.decompositionvar(a.getCoords(), c.getCoords());
			coeffC = c.decompositionvar(a.getCoords(), c.getCoords());
			coeffD = d.decompositionvar(a.getCoords(), c.getCoords());

			coeffAderiv = a.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau);
			coeffBderiv = b.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau);
			coeffCderiv = c.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau);
			coeffDderiv = d.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau);
		}

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
									ERI.LocalTwoCenterERIgd(A[i], B[j], C[k], D[l], tau);
						}

						if (coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] *
									ERI.LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}

						if (coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] *
									ERI.LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}

						if (coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] *
									ERI.LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}

						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] *
									ERI.LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}
					}
				}
			}
		}

		return sum * Constants.eV;
	}

	@Override
	public double Gg2d(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau1, int tau2) {
		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double[] coeffAderiv1 = a.derivativeDecomposition(a.getCoords(), c.getCoords(), tau1);
		double[] coeffBderiv1 = b.derivativeDecomposition(a.getCoords(), c.getCoords(), tau1);
		double[] coeffCderiv1 = c.derivativeDecomposition(a.getCoords(), c.getCoords(), tau1);
		double[] coeffDderiv1 = d.derivativeDecomposition(a.getCoords(), c.getCoords(), tau1);

		double[] coeffAderiv2 = a.derivativeDecomposition(a.getCoords(), c.getCoords(), tau2);
		double[] coeffBderiv2 = b.derivativeDecomposition(a.getCoords(), c.getCoords(), tau2);
		double[] coeffCderiv2 = c.derivativeDecomposition(a.getCoords(), c.getCoords(), tau2);
		double[] coeffDderiv2 = d.derivativeDecomposition(a.getCoords(), c.getCoords(), tau2);

		double[] coeffAderiv = a.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), tau1, tau2);
		double[] coeffBderiv = b.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), tau1, tau2);
		double[] coeffCderiv = c.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), tau1, tau2);
		double[] coeffDderiv = d.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), tau1, tau2);

		if (Math.abs(a.getCoords()[0] - c.getCoords()[0]) < 1E-3 &&
				Math.abs(a.getCoords()[1] - c.getCoords()[1]) < 1E-3) {

			coeffAderiv = a.secondDerivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1, tau2);
			coeffBderiv = b.secondDerivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1, tau2);
			coeffCderiv = c.secondDerivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1, tau2);
			coeffDderiv = d.secondDerivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1, tau2);

			coeffAderiv2 = a.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau2);
			coeffBderiv2 = b.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau2);
			coeffCderiv2 = c.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau2);
			coeffDderiv2 = d.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau2);

			coeffAderiv1 = a.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1);
			coeffBderiv1 = b.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1);
			coeffCderiv1 = c.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1);
			coeffDderiv1 = d.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1);

			coeffA = a.decompositionvar(a.getCoords(), c.getCoords());
			coeffB = b.decompositionvar(a.getCoords(), c.getCoords());
			coeffC = c.decompositionvar(a.getCoords(), c.getCoords());
			coeffD = d.decompositionvar(a.getCoords(), c.getCoords());
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
						double eri = ERI.LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						double erideriva = ERI.LocalTwoCenterERIgd(A[i], B[j], C[k], D[l], tau1);
						double eriderivb = ERI.LocalTwoCenterERIgd(A[i], B[j], C[k], D[l], tau2);
						double erideriv2 = ERI.LocalTwoCenterERIg2d(A[i], B[j], C[k], D[l], tau1, tau2);

						sum += coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] * eri;
						sum += coeffAderiv1[i] * coeffBderiv2[j] * coeffC[k] * coeffD[l] * eri;
						sum += coeffAderiv1[i] * coeffB[j] * coeffCderiv2[k] * coeffD[l] * eri;
						sum += coeffAderiv1[i] * coeffB[j] * coeffC[k] * coeffDderiv2[l] * eri;
						sum += coeffAderiv1[i] * coeffB[j] * coeffC[k] * coeffD[l] * eriderivb;

						sum += coeffAderiv2[i] * coeffBderiv1[j] * coeffC[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffBderiv1[j] * coeffCderiv2[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffBderiv1[j] * coeffC[k] * coeffDderiv2[l] * eri;
						sum += coeffA[i] * coeffBderiv1[j] * coeffC[k] * coeffD[l] * eriderivb;

						sum += coeffAderiv2[i] * coeffB[j] * coeffCderiv1[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffBderiv2[j] * coeffCderiv1[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffCderiv1[k] * coeffDderiv2[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffCderiv1[k] * coeffD[l] * eriderivb;

						sum += coeffAderiv2[i] * coeffB[j] * coeffC[k] * coeffDderiv1[l] * eri;
						sum += coeffA[i] * coeffBderiv2[j] * coeffC[k] * coeffDderiv1[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffCderiv2[k] * coeffDderiv1[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv1[l] * eriderivb;

						sum += coeffAderiv2[i] * coeffB[j] * coeffC[k] * coeffD[l] * erideriva;
						sum += coeffA[i] * coeffBderiv2[j] * coeffC[k] * coeffD[l] * erideriva;
						sum += coeffA[i] * coeffB[j] * coeffCderiv2[k] * coeffD[l] * erideriva;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv2[l] * erideriva;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] * erideriv2;
					}
				}
			}
		}

		return sum * Constants.eV;
	}

	@Override
	public double Gpd(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int num, int type) {
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
			p1deriv = a.getAtom().p1pd(type);
			p2deriv = a.getAtom().p2pd(type);
			D1deriv = a.getAtom().D1pd(type);
			D2deriv = a.getAtom().D2pd(type);
		}
		else if (num == 1) {
			p1deriv = c.getAtom().p1pd(type);
			p2deriv = c.getAtom().p2pd(type);
			D1deriv = c.getAtom().D1pd(type);
			D2deriv = c.getAtom().D2pd(type);
		}

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) sum2 +=
								coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
										ERI.LocalTwoCenterERIpd(A[i], B[j], C[k], D[l], D1deriv, D2deriv, p1deriv,
												p2deriv, num, type);
					}
				}
			}
		}

		return sum2 * Constants.eV;
	}

	@Override
	public double H(NDDO6G a, NDDO6G b) {
		return 0.5 * (a.beta + b.beta) * LCGTO.S(a, b);
	}

	@Override
	public double Hgd(NDDO6G a, NDDO6G b, int tau) {
		return 0.5 * (a.beta + b.beta) * LCGTO.Sgd(a, b, tau);
	}

	@Override
	public double Hg2d(NDDO6G a, NDDO6G b, int tau1, int tau2) {
		return 0.5 * (a.beta + b.beta) * LCGTO.Sg2d(a, b, tau1, tau2);
	}

	@Override
	public double Hzetapd(NDDO6G a, NDDO6G b, int num, int type) { // todo rename to betaparamzetaderiv
		if (num == -1) {
			return 0;
		}

		int hasA = a.getL() == type && num != 1 ? 1 : 0;
		int hasB = b.getL() == type && num != 0 ? 1 : 0;

		hasB <<= 1;
		int alltogether = hasA + hasB - 1;

		return 0.5 * (a.beta + b.beta) * STO6G.Spd(a, b, alltogether);
	}

	@Override
	public double Hbetapd(NDDO6G a, NDDO6G b, int num) {
		switch (num) {
			case 0:
			case 1:
				return 0.5 * LCGTO.S(a, b);
			case 2:
				return LCGTO.S(a, b);
			default:
				return 0;
		}
	}
}
