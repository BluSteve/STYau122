package nddo.defaults;

import nddo.Constants;
import nddo.NDDOOrbitalMethods;
import nddo.math.ERI;
import nddo.scf.LCGTO;
import nddo.scf.STO6G;

import static nddo.Constants.eV;
import static nddo.math.ERI.*;

public class NDDO6GMethods implements NDDOOrbitalMethods<NDDO6G> {
	private static double Gcrossp2d(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int num1, int type1, int num2,
									int type2) {
		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double sum2 = 0;

		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		double p11deriv = 0;
		double p21deriv = 0;
		double D11deriv = 0;
		double D21deriv = 0;

		double p12deriv = 0;
		double p22deriv = 0;
		double D12deriv = 0;
		double D22deriv = 0;

		if (num1 == 0 && num2 == 1) {
			p11deriv = a.getAtom().p1pd(type1);
			p21deriv = a.getAtom().p2pd(type1);
			D11deriv = a.getAtom().D1pd(type1);
			D21deriv = a.getAtom().D2pd(type1);

			p12deriv = c.getAtom().p1pd(type2);
			p22deriv = c.getAtom().p2pd(type2);
			D12deriv = c.getAtom().D1pd(type2);
			D22deriv = c.getAtom().D2pd(type2);
		}

		else if (num1 == 1 && num2 == 0) {
			p11deriv = a.getAtom().p1pd(type2);
			p21deriv = a.getAtom().p2pd(type2);
			D11deriv = a.getAtom().D1pd(type2);
			D21deriv = a.getAtom().D2pd(type2);

			p12deriv = c.getAtom().p1pd(type1);
			p22deriv = c.getAtom().p2pd(type1);
			D12deriv = c.getAtom().D1pd(type1);
			D22deriv = c.getAtom().D2pd(type1);
		}


		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {

						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum2 += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
									LocalTwoCenterERIcrossp2d(A[i], B[j], C[k], D[l], D11deriv, D21deriv, p11deriv,
											p21deriv, D12deriv, D22deriv, p12deriv, p22deriv) * eV;
						}
					}
				}
			}
		}


		return sum2;
	}

	private static double Gdiagp2d(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int num, int type1, int type2) {
		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double sum2 = 0;

		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		double p1deriva = 0;
		double p2deriva = 0;
		double D1deriva = 0;
		double D2deriva = 0;

		double p1derivb = 0;
		double p2derivb = 0;
		double D1derivb = 0;
		double D2derivb = 0;

		double p1deriv2 = 0;
		double p2deriv2 = 0;
		double D1deriv2 = 0;
		double D2deriv2 = 0;

		if (num == 0 || num == 2) {
			p1deriva = a.getAtom().p1pd(type1);
			p2deriva = a.getAtom().p2pd(type1);
			D1deriva = a.getAtom().D1pd(type1);
			D2deriva = a.getAtom().D2pd(type1);

			p1derivb = a.getAtom().p1pd(type2);
			p2derivb = a.getAtom().p2pd(type2);
			D1derivb = a.getAtom().D1pd(type2);
			D2derivb = a.getAtom().D2pd(type2);

			p1deriv2 = a.getAtom().p1p2d(type1 + type2);
			p2deriv2 = a.getAtom().p2p2d(type1 + type2);
			D1deriv2 = a.getAtom().D1p2d(type1 + type2);
			D2deriv2 = a.getAtom().D2p2d(type1 + type2);

		}
		else if (num == 1) {
			p1deriva = c.getAtom().p1pd(type1);
			p2deriva = c.getAtom().p2pd(type1);
			D1deriva = c.getAtom().D1pd(type1);
			D2deriva = c.getAtom().D2pd(type1);

			p1derivb = c.getAtom().p1pd(type2);
			p2derivb = c.getAtom().p2pd(type2);
			D1derivb = c.getAtom().D1pd(type2);
			D2derivb = c.getAtom().D2pd(type2);

			p1deriv2 = c.getAtom().p1p2d(type1 + type2);
			p2deriv2 = c.getAtom().p2p2d(type1 + type2);
			D1deriv2 = c.getAtom().D1p2d(type1 + type2);
			D2deriv2 = c.getAtom().D2p2d(type1 + type2);
		}

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
								sum2 += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
										LocalTwoCenterERIdiagp2d(A[i], B[j], C[k], D[l], D1deriva, D2deriva,
												p1deriva, p2deriva, D1derivb, D2derivb, p1derivb, p2derivb,
												D1deriv2, D2deriv2, p1deriv2, p2deriv2, num) * eV;
							}
						}
					}
				}
			}
		}

		return sum2;
	}

	private static double Gcrossp2gd(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int num1, int type1, int num2,
									int type2, int tau) {


		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		double p11deriv = 0;
		double p21deriv = 0;
		double D11deriv = 0;
		double D21deriv = 0;

		double p12deriv = 0;
		double p22deriv = 0;
		double D12deriv = 0;
		double D22deriv = 0;

		if (num1 == 0 && num2 == 1) {
			p11deriv = a.getAtom().p1pd(type1);
			p21deriv = a.getAtom().p2pd(type1);
			D11deriv = a.getAtom().D1pd(type1);
			D21deriv = a.getAtom().D2pd(type1);

			p12deriv = c.getAtom().p1pd(type2);
			p22deriv = c.getAtom().p2pd(type2);
			D12deriv = c.getAtom().D1pd(type2);
			D22deriv = c.getAtom().D2pd(type2);
		}

		else if (num1 == 1 && num2 == 0) {
			p11deriv = a.getAtom().p1pd(type2);
			p21deriv = a.getAtom().p2pd(type2);
			D11deriv = a.getAtom().D1pd(type2);
			D21deriv = a.getAtom().D2pd(type2);

			p12deriv = c.getAtom().p1pd(type1);
			p22deriv = c.getAtom().p2pd(type1);
			D12deriv = c.getAtom().D1pd(type1);
			D22deriv = c.getAtom().D2pd(type1);
		}

		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double[] coeffAderiv = a.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffBderiv = b.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffCderiv = c.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffDderiv = d.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);

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

		double sum = 0;

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) sum +=
								coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
										LocalTwoCenterERIcrossp2gd(A[i], B[j], C[k], D[l], D11deriv, D21deriv,
												p11deriv, p21deriv, D12deriv, D22deriv, p12deriv, p22deriv, tau);

						if (coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) sum +=
								coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] *
										LocalTwoCenterERIcrossp2d(A[i], B[j], C[k], D[l], D11deriv, D21deriv, p11deriv,

												p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
						if (coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] != 0) sum +=
								coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] *
										LocalTwoCenterERIcrossp2d(A[i], B[j], C[k], D[l], D11deriv, D21deriv, p11deriv,
												p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);

						if (coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] != 0) sum +=
								coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] *
										LocalTwoCenterERIcrossp2d(A[i], B[j], C[k], D[l], D11deriv, D21deriv, p11deriv,
												p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);

						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] != 0) sum +=
								coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] *
										LocalTwoCenterERIcrossp2d(A[i], B[j], C[k], D[l], D11deriv, D21deriv, p11deriv,
												p21deriv, D12deriv, D22deriv, p12deriv, p22deriv);
					}
				}
			}
		}

		return sum * Constants.eV;
	}

	private static double Gdiagp2gd(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int num, int type1, int type2,
								   int tau) {
		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		double p1deriva = 0;
		double p2deriva = 0;
		double D1deriva = 0;
		double D2deriva = 0;

		double p1derivb = 0;
		double p2derivb = 0;
		double D1derivb = 0;
		double D2derivb = 0;

		double p1deriv2 = 0;
		double p2deriv2 = 0;
		double D1deriv2 = 0;
		double D2deriv2 = 0;

		if (num == 0 || num == 2) {
			p1deriva = a.getAtom().p1pd(type1);
			p2deriva = a.getAtom().p2pd(type1);
			D1deriva = a.getAtom().D1pd(type1);
			D2deriva = a.getAtom().D2pd(type1);

			p1derivb = a.getAtom().p1pd(type2);
			p2derivb = a.getAtom().p2pd(type2);
			D1derivb = a.getAtom().D1pd(type2);
			D2derivb = a.getAtom().D2pd(type2);

			p1deriv2 = a.getAtom().p1p2d(type1 + type2);
			p2deriv2 = a.getAtom().p2p2d(type1 + type2);
			D1deriv2 = a.getAtom().D1p2d(type1 + type2);
			D2deriv2 = a.getAtom().D2p2d(type1 + type2);

		}
		else if (num == 1) {
			p1deriva = c.getAtom().p1pd(type1);
			p2deriva = c.getAtom().p2pd(type1);
			D1deriva = c.getAtom().D1pd(type1);
			D2deriva = c.getAtom().D2pd(type1);

			p1derivb = c.getAtom().p1pd(type2);
			p2derivb = c.getAtom().p2pd(type2);
			D1derivb = c.getAtom().D1pd(type2);
			D2derivb = c.getAtom().D2pd(type2);

			p1deriv2 = c.getAtom().p1p2d(type1 + type2);
			p2deriv2 = c.getAtom().p2p2d(type1 + type2);
			D1deriv2 = c.getAtom().D1p2d(type1 + type2);
			D2deriv2 = c.getAtom().D2p2d(type1 + type2);
		}

		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double[] coeffAderiv = a.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffBderiv = b.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffCderiv = c.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffDderiv = d.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);

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
		double sum = 0;

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
									LocalTwoCenterERIdiagp2gd(A[i], B[j], C[k], D[l], D1deriva, D2deriva, p1deriva,
											p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2,
											p1deriv2, p2deriv2, num, tau);
						}

						double erideriv2 = LocalTwoCenterERIdiagp2d(A[i], B[j], C[k], D[l], D1deriva, D2deriva,
								p1deriva, p2deriva, D1derivb, D2derivb, p1derivb, p2derivb, D1deriv2, D2deriv2,
								p1deriv2, p2deriv2, num);

						if (coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] * erideriv2;
						}

						if (coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] * erideriv2;
						}

						if (coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] * erideriv2;
						}

						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] * erideriv2;
						}
					}
				}
			}
		}


		return sum * Constants.eV;
	}

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
									LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}

		return sum2 * eV;
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
									LocalTwoCenterERIgd(A[i], B[j], C[k], D[l], tau);
						}

						if (coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] *
									LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}

						if (coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] *
									LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}

						if (coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] *
									LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}

						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] *
									LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						}
					}
				}
			}
		}

		return sum * eV;
	}

	@Override
	public double Gg2d(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau1, int tau2) {
		double[] coeffA;
		double[] coeffB;
		double[] coeffC;
		double[] coeffD;

		double[] coeffAderiva;
		double[] coeffBderiva;
		double[] coeffCderiva;
		double[] coeffDderiva;

		double[] coeffAderivb;
		double[] coeffBderivb;
		double[] coeffCderivb;
		double[] coeffDderivb;

		double[] coeffAderiv2;
		double[] coeffBderiv2;
		double[] coeffCderiv2;
		double[] coeffDderiv2;

		if (Math.abs(a.getCoords()[0] - c.getCoords()[0]) < 1E-3 &&
				Math.abs(a.getCoords()[1] - c.getCoords()[1]) < 1E-3) {

			coeffAderiv2 = a.secondDerivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1, tau2);
			coeffBderiv2 = b.secondDerivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1, tau2);
			coeffCderiv2 = c.secondDerivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1, tau2);
			coeffDderiv2 = d.secondDerivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1, tau2);

			coeffAderivb = a.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau2);
			coeffBderivb = b.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau2);
			coeffCderivb = c.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau2);
			coeffDderivb = d.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau2);

			coeffAderiva = a.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1);
			coeffBderiva = b.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1);
			coeffCderiva = c.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1);
			coeffDderiva = d.derivativeDecompositionvar(a.getCoords(), c.getCoords(), tau1);

			coeffA = a.decompositionvar(a.getCoords(), c.getCoords());
			coeffB = b.decompositionvar(a.getCoords(), c.getCoords());
			coeffC = c.decompositionvar(a.getCoords(), c.getCoords());
			coeffD = d.decompositionvar(a.getCoords(), c.getCoords());
		}

		else {
			coeffA = a.decomposition(a.getCoords(), c.getCoords());
			coeffB = b.decomposition(a.getCoords(), c.getCoords());
			coeffC = c.decomposition(a.getCoords(), c.getCoords());
			coeffD = d.decomposition(a.getCoords(), c.getCoords());

			coeffAderiva = a.derivativeDecomposition(a.getCoords(), c.getCoords(), tau1);
			coeffBderiva = b.derivativeDecomposition(a.getCoords(), c.getCoords(), tau1);
			coeffCderiva = c.derivativeDecomposition(a.getCoords(), c.getCoords(), tau1);
			coeffDderiva = d.derivativeDecomposition(a.getCoords(), c.getCoords(), tau1);

			coeffAderivb = a.derivativeDecomposition(a.getCoords(), c.getCoords(), tau2);
			coeffBderivb = b.derivativeDecomposition(a.getCoords(), c.getCoords(), tau2);
			coeffCderivb = c.derivativeDecomposition(a.getCoords(), c.getCoords(), tau2);
			coeffDderivb = d.derivativeDecomposition(a.getCoords(), c.getCoords(), tau2);

			coeffAderiv2 = a.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), tau1, tau2);
			coeffBderiv2 = b.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), tau1, tau2);
			coeffCderiv2 = c.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), tau1, tau2);
			coeffDderiv2 = d.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), tau1, tau2);
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

						double eri = LocalTwoCenterERI(A[i], B[j], C[k], D[l]);
						double erideriva = LocalTwoCenterERIgd(A[i], B[j], C[k], D[l], tau1);
						double eriderivb = LocalTwoCenterERIgd(A[i], B[j], C[k], D[l], tau2);
						double erideriv2 = LocalTwoCenterERIg2d(A[i], B[j], C[k], D[l], tau1, tau2);

						sum += coeffAderiv2[i] * coeffB[j] * coeffC[k] * coeffD[l] * eri;
						sum += coeffAderiva[i] * coeffBderivb[j] * coeffC[k] * coeffD[l] * eri;
						sum += coeffAderiva[i] * coeffB[j] * coeffCderivb[k] * coeffD[l] * eri;
						sum += coeffAderiva[i] * coeffB[j] * coeffC[k] * coeffDderivb[l] * eri;
						sum += coeffAderiva[i] * coeffB[j] * coeffC[k] * coeffD[l] * eriderivb;

						sum += coeffAderivb[i] * coeffBderiva[j] * coeffC[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffBderiv2[j] * coeffC[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffBderiva[j] * coeffCderivb[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffBderiva[j] * coeffC[k] * coeffDderivb[l] * eri;
						sum += coeffA[i] * coeffBderiva[j] * coeffC[k] * coeffD[l] * eriderivb;

						sum += coeffAderivb[i] * coeffB[j] * coeffCderiva[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffBderivb[j] * coeffCderiva[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffCderiv2[k] * coeffD[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffCderiva[k] * coeffDderivb[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffCderiva[k] * coeffD[l] * eriderivb;

						sum += coeffAderivb[i] * coeffB[j] * coeffC[k] * coeffDderiva[l] * eri;
						sum += coeffA[i] * coeffBderivb[j] * coeffC[k] * coeffDderiva[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffCderivb[k] * coeffDderiva[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv2[l] * eri;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiva[l] * eriderivb;

						sum += coeffAderivb[i] * coeffB[j] * coeffC[k] * coeffD[l] * erideriva;
						sum += coeffA[i] * coeffBderivb[j] * coeffC[k] * coeffD[l] * erideriva;
						sum += coeffA[i] * coeffB[j] * coeffCderivb[k] * coeffD[l] * erideriva;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderivb[l] * erideriva;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] * erideriv2;
					}
				}
			}
		}

		return sum * eV;
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
										LocalTwoCenterERIpd(A[i], B[j], C[k], D[l], D1deriv, D2deriv, p1deriv,
												p2deriv, num);
					}
				}
			}
		}

		return sum2 * eV;
	}

	@Override
	public double Gp2d(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int num1, int type1, int num2, int type2) {
		if (num1 == -1 || num2 == -1) {
			return 0;
		}

		else if (num1 + num2 == 1) {
			return Gcrossp2d(a, b, c, d, num1, type1, num2, type2);
		}
		else if (num1 == num2) {
			return Gdiagp2d(a, b, c, d, num1, type1, type2);
		}

		return 0;
	}

	@Override
	public double Gpgd(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int num, int type, int tau) {
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

		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double[] coeffAderiv = a.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffBderiv = b.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffCderiv = c.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);
		double[] coeffDderiv = d.derivativeDecomposition(a.getCoords(), c.getCoords(), tau);

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

		double sum = 0;

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
									LocalTwoCenterERIpgd(A[i], B[j], C[k], D[l], D1deriv, D2deriv, p1deriv, p2deriv,
											num, tau);
						}

						double erideriv = LocalTwoCenterERIpd(A[i], B[j], C[k], D[l], D1deriv,
								D2deriv, p1deriv, p2deriv, num);

						if (coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] * erideriv;
						}

						if (coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] * erideriv;
						}

						if (coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] * erideriv;
						}

						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] * erideriv;
						}
					}
				}
			}
		}


		return sum * eV;
	}

	@Override
	public double Gp2gd(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int num1,  int type1, int num2, int type2, int tau) {
		if (num1 == -1 || num2 == -1) {
			return 0;
		}

		else if (num1 + num2 == 1) {
			return Gcrossp2gd(a, b, c, d, num1, type1, num2, type2, tau);
		}
		else if (num1 == num2) {
			return Gdiagp2gd(a, b, c, d, num1, type1, type2, tau);
		}

		return 0;
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
	public double Hzetapd(NDDO6G a, NDDO6G b, int num, int type) {
		return 0.5 * (a.beta + b.beta) * STO6G.Spd(a, b, num, type);
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

	@Override
	public double Hzetazetap2d(NDDO6G a, NDDO6G b, int num1, int type1, int num2, int type2) {
		return 0.5 * (a.beta + b.beta) * STO6G.Sp2d(a, b, num1, type1, num2, type2);
	}

	@Override
	public double Hbetazetap2d(NDDO6G a, NDDO6G b, int num1, int num2, int type2) {
		switch (num1) {
			case 0:
			case 1:
				return 0.5 * STO6G.Spd(a, b, num2, type2);
			case 2:
				return STO6G.Spd(a, b, num2, type2);
			default:
				return 0;
		}
	}

	@Override
	public double Hzetapgd(NDDO6G a, NDDO6G b, int num, int type, int tau) {
		return 0.5 * (a.beta + b.beta) * STO6G.Spgd(a, b, num, type, tau);
	}

	@Override
	public double Hbetapgd(NDDO6G a, NDDO6G b, int num, int tau) {
		switch (num) {
			case 0:
			case 1:
				return 0.5 * STO6G.Sgd(a, b, tau);
			case 2:
				return STO6G.Sgd(a, b, tau);
			default:
				return 0;
		}
	}

	@Override
	public double Hzetazetap2gd(NDDO6G a, NDDO6G b, int num1, int type1, int num2, int type2, int tau) {
		return 0.5 * (a.beta + b.beta) * STO6G.Sp2gd(a, b, num1, type1, num2, type2, tau);
	}

	@Override
	public double Hbetazetap2gd(NDDO6G a, NDDO6G b, int num1, int num2, int type2, int tau) {
		switch (num1) {
			case 0:
			case 1:
				return 0.5 * STO6G.Spgd(a, b, num2, type2, tau);
			case 2:
				return STO6G.Spgd(a, b, num2, type2, tau);
			default:
				return 0;
		}
	}
}
