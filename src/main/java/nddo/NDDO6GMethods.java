package nddo;

import nddo.math.ERI;
import nddo.scf.LCGTO;
import nddo.scf.Orbital;
import nddo.scf.STO6G;

public class NDDO6GMethods implements NDDOOrbitalMethods<NDDO6G> {

	@Override
	public double OneCenterERI(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d) {
		return ERI.OneCenterERI(a,b,c,d);
	}

	@Override
	public double getG(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d) {
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
									ERI.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * Constants.eV;
						}

		return sum2;
	}

	@Override
	public double getGderiv(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d,
							int tau) {


		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		if (Math.abs(a.getCoords()[0] - c.getCoords()[0]) < 1E-3 &&
				Math.abs(a.getCoords()[1] - c.getCoords()[1]) < 1E-3) {

			double[] coeffA2 = a.decomposition2(a.getCoords(), c.getCoords());
			double[] coeffB2 = b.decomposition2(a.getCoords(), c.getCoords());
			double[] coeffC2 = c.decomposition2(a.getCoords(), c.getCoords());
			double[] coeffD2 = d.decomposition2(a.getCoords(), c.getCoords());

			double[] coeffAderiv2 = Orbital.derivativeDecomposition2(a.getCoords(), c.getCoords(), a, tau);
			double[] coeffBderiv2 = Orbital.derivativeDecomposition2(a.getCoords(), c.getCoords(), b, tau);
			double[] coeffCderiv2 = Orbital.derivativeDecomposition2(a.getCoords(), c.getCoords(), c, tau);
			double[] coeffDderiv2 = Orbital.derivativeDecomposition2(a.getCoords(), c.getCoords(), d, tau);


			double sum2 = 0;

			for (int i = 0; i < coeffA2.length; i++) {
				for (int j = 0; j < coeffB2.length; j++) {
					for (int k = 0; k < coeffC2.length; k++) {
						for (int l = 0; l < coeffD2.length; l++) {


							if (coeffA2[i] * coeffB2[j] * coeffC2[k] *
									coeffD2[l] != 0) {
								sum2 += coeffA2[i] * coeffB2[j] * coeffC2[k] *
										coeffD2[l] *
										ERI.LocalTwoCenterERIderiv(A[i], B[j],
												C[k],
												D[l], tau) * Constants.eV;
							}
							if (coeffAderiv2[i] * coeffB2[j] * coeffC2[k] *
									coeffD2[l] !=
									0) {
								sum2 += coeffAderiv2[i] * coeffB2[j] *
										coeffC2[k] *
										coeffD2[l] *
										ERI.LocalTwoCenterERI(A[i], B[j],
												C[k], D[l]) *
										Constants.eV;
							}

							if (coeffA2[i] * coeffBderiv2[j] * coeffC2[k] *
									coeffD2[l] !=
									0) {
								sum2 += coeffA2[i] * coeffBderiv2[j] *
										coeffC2[k] *
										coeffD2[l] *
										ERI.LocalTwoCenterERI(A[i], B[j],
												C[k], D[l]) *
										Constants.eV;
							}

							if (coeffA2[i] * coeffB2[j] * coeffCderiv2[k] *
									coeffD2[l] !=
									0) {
								sum2 += coeffA2[i] * coeffB2[j] *
										coeffCderiv2[k] *
										coeffD2[l] *
										ERI.LocalTwoCenterERI(A[i], B[j],
												C[k], D[l]) *
										Constants.eV;
							}

							if (coeffA2[i] * coeffB2[j] * coeffC2[k] *
									coeffDderiv2[l] !=
									0) {
								sum2 += coeffA2[i] * coeffB2[j] * coeffC2[k] *
										coeffDderiv2[l] *
										ERI.LocalTwoCenterERI(A[i], B[j],
												C[k], D[l]) *
										Constants.eV;
							}


						}
					}
				}
			}

			return sum2;
		}
		else {
			double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
			double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
			double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
			double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

			double[] coeffAderiv =
					Orbital.derivativeDecomposition(a.getCoords(), c.getCoords(), a,
							tau);
			double[] coeffBderiv =
					Orbital.derivativeDecomposition(a.getCoords(), c.getCoords(), b,
							tau);
			double[] coeffCderiv =
					Orbital.derivativeDecomposition(a.getCoords(), c.getCoords(), c,
							tau);
			double[] coeffDderiv =
					Orbital.derivativeDecomposition(a.getCoords(), c.getCoords(), d,
							tau);


			double sum = 0;

			for (int i = 0; i < coeffA.length; i++) {
				for (int j = 0; j < coeffB.length; j++) {
					for (int k = 0; k < coeffC.length; k++) {
						for (int l = 0; l < coeffD.length; l++) {


							if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] !=
									0) {
								sum += coeffA[i] * coeffB[j] * coeffC[k] *
										coeffD[l] *
										ERI.LocalTwoCenterERIderiv(A[i], B[j],
												C[k],
												D[l], tau) * Constants.eV;
							}
							if (coeffAderiv[i] * coeffB[j] * coeffC[k] *
									coeffD[l] != 0) {
								sum += coeffAderiv[i] * coeffB[j] * coeffC[k] *
										coeffD[l] *
										ERI.LocalTwoCenterERI(A[i], B[j],
												C[k], D[l]) *
										Constants.eV;
							}

							if (coeffA[i] * coeffBderiv[j] * coeffC[k] *
									coeffD[l] != 0) {
								sum += coeffA[i] * coeffBderiv[j] * coeffC[k] *
										coeffD[l] *
										ERI.LocalTwoCenterERI(A[i], B[j],
												C[k], D[l]) *
										Constants.eV;
							}

							if (coeffA[i] * coeffB[j] * coeffCderiv[k] *
									coeffD[l] != 0) {
								sum += coeffA[i] * coeffB[j] * coeffCderiv[k] *
										coeffD[l] *
										ERI.LocalTwoCenterERI(A[i], B[j],
												C[k], D[l]) *
										Constants.eV;
							}

							if (coeffA[i] * coeffB[j] * coeffC[k] *
									coeffDderiv[l] != 0) {
								sum += coeffA[i] * coeffB[j] * coeffC[k] *
										coeffDderiv[l] *
										ERI.LocalTwoCenterERI(A[i], B[j],
												C[k], D[l]) *
										Constants.eV;
							}


						}
					}
				}
			}


			return sum;
		}
	}

	@Override
	public double getGderiv2(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d,
							 int tau1, int tau2) {


		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

		double[] coeffAderiv1 = Orbital
				.derivativeDecomposition(a.getCoords(), c.getCoords(), a,
						tau1);
		double[] coeffBderiv1 = Orbital
				.derivativeDecomposition(a.getCoords(), c.getCoords(), b,
						tau1);
		double[] coeffCderiv1 = Orbital
				.derivativeDecomposition(a.getCoords(), c.getCoords(), c,
						tau1);
		double[] coeffDderiv1 = Orbital
				.derivativeDecomposition(a.getCoords(), c.getCoords(), d,
						tau1);

		double[] coeffAderiv2 = Orbital
				.derivativeDecomposition(a.getCoords(), c.getCoords(), a,
						tau2);
		double[] coeffBderiv2 = Orbital
				.derivativeDecomposition(a.getCoords(), c.getCoords(), b,
						tau2);
		double[] coeffCderiv2 = Orbital
				.derivativeDecomposition(a.getCoords(), c.getCoords(), c,
						tau2);
		double[] coeffDderiv2 = Orbital
				.derivativeDecomposition(a.getCoords(), c.getCoords(), d,
						tau2);

		double[] coeffAderiv =
				Orbital.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), a,
						tau1, tau2);
		double[] coeffBderiv =
				Orbital.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), b,
						tau1, tau2);
		double[] coeffCderiv =
				Orbital.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), c,
						tau1, tau2);
		double[] coeffDderiv =
				Orbital.secondDerivativeDecomposition(a.getCoords(), c.getCoords(), d,
						tau1, tau2);

		if (Math.abs(a.getCoords()[0] - c.getCoords()[0]) < 1E-3 &&
				Math.abs(a.getCoords()[1] - c.getCoords()[1]) < 1E-3) {

			coeffAderiv =
					Orbital.secondDerivativeDecomposition2(a.getCoords(),
							c.getCoords(),
							a, tau1, tau2);
			coeffBderiv =
					Orbital.secondDerivativeDecomposition2(a.getCoords(),
							c.getCoords(),
							b, tau1, tau2);
			coeffCderiv =
					Orbital.secondDerivativeDecomposition2(a.getCoords(),
							c.getCoords(),
							c, tau1, tau2);
			coeffDderiv =
					Orbital.secondDerivativeDecomposition2(a.getCoords(),
							c.getCoords(),
							d, tau1, tau2);

			coeffAderiv2 = Orbital
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), a,
							tau2);
			coeffBderiv2 = Orbital
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), b,
							tau2);
			coeffCderiv2 = Orbital
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), c,
							tau2);
			coeffDderiv2 = Orbital
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), d,
							tau2);

			coeffAderiv1 = Orbital
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), a,
							tau1);
			coeffBderiv1 = Orbital
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), b,
							tau1);
			coeffCderiv1 = Orbital
					.derivativeDecomposition2(a.getCoords(), c.getCoords(), c,
							tau1);
			coeffDderiv1 = Orbital
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

						double eri = ERI.LocalTwoCenterERI(A[i], B[j], C[k],
								D[l]);

						double erideriv1 = ERI
								.LocalTwoCenterERIderiv(A[i], B[j], C[k], D[l],
										tau1);

						double erideriv2 = ERI
								.LocalTwoCenterERIderiv(A[i], B[j], C[k], D[l],
										tau2);

						double erideriv =
								ERI.LocalTwoCenterERIderiv2(A[i], B[j], C[k], D[l],
										tau1, tau2);

						sum += coeffAderiv[i] * coeffB[j] * coeffC[k] *
								coeffD[l] * eri * Constants.eV;
						sum += coeffAderiv1[i] * coeffBderiv2[j] * coeffC[k] *
								coeffD[l] * eri * Constants.eV;
						sum += coeffAderiv1[i] * coeffB[j] * coeffCderiv2[k] *
								coeffD[l] * eri * Constants.eV;
						sum += coeffAderiv1[i] * coeffB[j] * coeffC[k] *
								coeffDderiv2[l] * eri * Constants.eV;
						sum += coeffAderiv1[i] * coeffB[j] * coeffC[k] *
								coeffD[l] * erideriv2 * Constants.eV;

						sum += coeffAderiv2[i] * coeffBderiv1[j] * coeffC[k] *
								coeffD[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffBderiv[j] * coeffC[k] *
								coeffD[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffBderiv1[j] * coeffCderiv2[k] *
								coeffD[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffBderiv1[j] * coeffC[k] *
								coeffDderiv2[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffBderiv1[j] * coeffC[k] *
								coeffD[l] * erideriv2 * Constants.eV;

						sum += coeffAderiv2[i] * coeffB[j] * coeffCderiv1[k] *
								coeffD[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffBderiv2[j] * coeffCderiv1[k] *
								coeffD[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffB[j] * coeffCderiv[k] *
								coeffD[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffB[j] * coeffCderiv1[k] *
								coeffDderiv2[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffB[j] * coeffCderiv1[k] *
								coeffD[l] * erideriv2 * Constants.eV;

						sum += coeffAderiv2[i] * coeffB[j] * coeffC[k] *
								coeffDderiv1[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffBderiv2[j] * coeffC[k] *
								coeffDderiv1[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffB[j] * coeffCderiv2[k] *
								coeffDderiv1[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffB[j] * coeffC[k] *
								coeffDderiv[l] * eri * Constants.eV;
						sum += coeffA[i] * coeffB[j] * coeffC[k] *
								coeffDderiv1[l] * erideriv2 * Constants.eV;

						sum += coeffAderiv2[i] * coeffB[j] * coeffC[k] *
								coeffD[l] * erideriv1 * Constants.eV;
						sum += coeffA[i] * coeffBderiv2[j] * coeffC[k] *
								coeffD[l] * erideriv1 * Constants.eV;
						sum += coeffA[i] * coeffB[j] * coeffCderiv2[k] *
								coeffD[l] * erideriv1 * Constants.eV;
						sum += coeffA[i] * coeffB[j] * coeffC[k] *
								coeffDderiv2[l] * erideriv1 * Constants.eV;
						sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
								erideriv * Constants.eV;


					}
				}
			}
		}

		return sum;
	}

	@Override
	public double getGderiv(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d,
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
			p1deriv = a.getAtom().p1Deriv(type);
			p2deriv = a.getAtom().p2Deriv(type);
			D1deriv = a.getAtom().D1Deriv(type);
			D2deriv = a.getAtom().D2Deriv(type);
		}
		else if (num == 1) {
			p1deriv = c.getAtom().p1Deriv(type);
			p2deriv = c.getAtom().p2Deriv(type);
			D1deriv = c.getAtom().D1Deriv(type);
			D2deriv = c.getAtom().D2Deriv(type);
		}

		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {

						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] !=
								0) {
							sum2 += coeffA[i] * coeffB[j] * coeffC[k] *
									coeffD[l] *
									ERI.LocalTwoCenterERIderiv(A[i], B[j], C[k],
											D[l], D1deriv, D2deriv, p1deriv,
											p2deriv, num, type) * Constants.eV;
						}


					}
				}
			}
		}


		return sum2;
	}

	@Override
	public double beta(NDDO6G a, NDDO6G b) {
		return 0.5 * (a.beta + b.beta) * LCGTO.getS(a, b);
	}

	@Override
	public double betaderiv(NDDO6G a, NDDO6G b, int tau) {
		return 0.5 * (a.beta + b.beta) * LCGTO.getSDeriv(a, b, tau);
	}

	@Override
	public double betaderiv2(NDDO6G a, NDDO6G b, int tau1, int tau2) {
		return 0.5 * (a.beta + b.beta) * LCGTO.getSDeriv2(a, b, tau1, tau2);
	}

	@Override
	public double betaparamderiv(NDDO6G a, NDDO6G b, int num, int type) { // todo rename to betaparamzetaderiv

		if (num == -1) {
			return 0;
		}

		int hasA = (a.getL() == type) && num != 1 ? 1 : 0;
		int hasB = (b.getL() == type) && num != 0 ? 1 : 0;

		hasB <<= 1;
		int alltogether = hasA + hasB - 1;

		return 0.5 * (a.beta + b.beta) * STO6G.getSderivzeta(a, b, alltogether);
	}

	@Override
	public double betaparambetaderiv(NDDO6G a, NDDO6G b, double sum) {
		return sum * LCGTO.getS(a, b);
	}
}
