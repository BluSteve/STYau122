package nddo.geometry;

import nddo.Constants;
import nddo.NDDO6G;
import nddo.NDDOAtom;
import nddo.math.ERI;
import nddo.scf.GTO;
import nddo.scf.LCGTO;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

public class GeometryDerivative {

	public static double getGderiv(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d,
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

			double[] coeffAderiv2 = derivativeDecomposition2(a.getCoords(), c.getCoords(), a, tau);
			double[] coeffBderiv2 = derivativeDecomposition2(a.getCoords(), c.getCoords(), b, tau);
			double[] coeffCderiv2 = derivativeDecomposition2(a.getCoords(), c.getCoords(), c, tau);
			double[] coeffDderiv2 = derivativeDecomposition2(a.getCoords(), c.getCoords(), d, tau);


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
					derivativeDecomposition(a.getCoords(), c.getCoords(), a,
							tau);
			double[] coeffBderiv =
					derivativeDecomposition(a.getCoords(), c.getCoords(), b,
							tau);
			double[] coeffCderiv =
					derivativeDecomposition(a.getCoords(), c.getCoords(), c,
							tau);
			double[] coeffDderiv =
					derivativeDecomposition(a.getCoords(), c.getCoords(), d,
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

	public static double getGderivfinite(NDDO6G a, NDDO6G b, NDDO6G c,
										 NDDO6G d,
										 int tau) {

		double orig = NDDO6G.getG(a, b, c, d);

		double[] newcoords = a.getCoords().clone();

		newcoords[tau] += 1E-8;

		NDDO6G anew = new NDDO6G(a, newcoords);
		NDDO6G bnew = new NDDO6G(b, newcoords);

		double perturbed = NDDO6G.getG(anew, bnew, c, d);

		return (perturbed - orig) / 1E-8;
	}


	public static double[] derivativeDecomposition(double[] point1,
												   double[] point2, NDDO6G a,
												   int tau) {

		if (a.getL() == 0) {
			return new double[]{0};
		}

		double R = GTO.R(point1, point2);
		double Rxy =
				Math.sqrt((point2[1] - point1[1]) * (point2[1] - point1[1]) +
						(point2[0] - point1[0]) * (point2[0] - point1[0]));


		switch (tau) {
			case 0:
				if (a.geti() == 1) {
					double x1 = (point2[2] - point1[2]) / (R * Rxy) -
							(point2[0] - point1[0]) * (point2[0] - point1[0]) *
									(point2[2] - point1[2]) /
									(Rxy * Rxy * Rxy * R) -
							(point2[0] - point1[0]) * (point2[0] - point1[0]) *
									(point2[2] - point1[2]) / (R * R * R * Rxy);
					double x2 =
							-(point2[0] - point1[0]) * (point2[1] - point1[1]) /
									(Rxy * Rxy * Rxy);
					double x3 =
							(point2[0] - point1[0]) * (point2[0] - point1[0]) /
									(R * R * R) - 1 / R;

					return new double[]{x1, x2, x3};
				}
				else if (a.getj() == 1) {
					double x1 =
							-(point2[0] - point1[0]) * (point2[1] - point1[1]) *
									(point2[2] - point1[2]) /
									(Rxy * Rxy * Rxy * R) -
									(point2[0] - point1[0]) *
											(point2[1] - point1[1]) *
											(point2[2] - point1[2]) /
											(R * R * R * Rxy);
					double x2 =
							(point2[0] - point1[0]) * (point2[0] - point1[0]) /
									(Rxy * Rxy * Rxy) - 1 / Rxy;
					double x3 =
							(point2[0] - point1[0]) * (point2[1] - point1[1]) /
									(R * R * R);

					return new double[]{x1, x2, x3};
				}
				else if (a.getk() == 1) {
					double x1 = (point2[0] - point1[0]) * Rxy / (R * R * R) -
							(point2[0] - point1[0]) / (R * Rxy);
					double x2 = 0;
					double x3 =
							(point2[0] - point1[0]) * (point2[2] - point1[2]) /
									(R * R * R);

					return new double[]{x1, x2, x3};
				}
			case 1:
				if (a.geti() == 1) {
					double x1 =
							-(point2[0] - point1[0]) * (point2[1] - point1[1]) *
									(point2[2] - point1[2]) /
									(Rxy * Rxy * Rxy * R) -
									(point2[0] - point1[0]) *
											(point2[1] - point1[1]) *
											(point2[2] - point1[2]) /
											(R * R * R * Rxy);
					double x2 =
							-(point2[1] - point1[1]) * (point2[1] - point1[1]) /
									(Rxy * Rxy * Rxy) + 1 / Rxy;
					double x3 =
							(point2[0] - point1[0]) * (point2[1] - point1[1]) /
									(R * R * R);

					return new double[]{x1, x2, x3};
				}
				else if (a.getj() == 1) {
					double x1 = (point2[2] - point1[2]) / (R * Rxy) -
							(point2[2] - point1[2]) * (point2[1] - point1[1]) *
									(point2[1] - point1[1]) /
									(Rxy * Rxy * Rxy * R) -
							(point2[2] - point1[2]) * (point2[1] - point1[1]) *
									(point2[1] - point1[1]) / (R * R * R * Rxy);
					double x2 =
							(point2[0] - point1[0]) * (point2[1] - point1[1]) /
									(Rxy * Rxy * Rxy);
					double x3 =
							(point2[1] - point1[1]) * (point2[1] - point1[1]) /
									(R * R * R) - 1 / R;

					return new double[]{x1, x2, x3};
				}
				else if (a.getk() == 1) {
					double x1 = (point2[1] - point1[1]) * Rxy / (R * R * R) -
							(point2[1] - point1[1]) / (R * Rxy);
					double x2 = 0;
					double x3 =
							(point2[2] - point1[2]) * (point2[1] - point1[1]) /
									(R * R * R);

					return new double[]{x1, x2, x3};
				}
			case 2:
				if (a.geti() == 1) {
					double x1 = (point2[0] - point1[0]) / (R * Rxy) -
							(point2[0] - point1[0]) * (point2[2] - point1[2]) *
									(point2[2] - point1[2]) / (R * R * R * Rxy);
					double x2 = 0;
					double x3 =
							(point2[2] - point1[2]) * (point2[0] - point1[0]) /
									(R * R * R);

					return new double[]{x1, x2, x3};
				}
				else if (a.getj() == 1) {
					double x1 = (point2[1] - point1[1]) / (R * Rxy) -
							(point2[1] - point1[1]) * (point2[2] - point1[2]) *
									(point2[2] - point1[2]) / (R * R * R * Rxy);
					double x2 = 0;
					double x3 =
							(point2[2] - point1[2]) * (point2[1] - point1[1]) /
									(R * R * R);

					return new double[]{x1, x2, x3};
				}
				else if (a.getk() == 1) {
					double x1 = (point2[2] - point1[2]) * Rxy / (R * R * R);
					double x2 = 0;
					double x3 =
							(point2[2] - point1[2]) * (point2[2] - point1[2]) /
									(R * R * R) - 1 / R;

					return new double[]{x1, x2, x3};
				}
		}
		return null;
	}

	public static double[] derivativeDecomposition2(double[] point1,
													double[] point2, NDDO6G a,
													int tau) {

		if (a.getL() == 0) {
			return new double[]{0};
		}

		double x = point2[0] - point1[0];

		double y = point2[1] - point1[1];

		double z = point2[2] - point1[2];

		double Rxz = Math.sqrt(x * x + z * z);

		double R = GTO.R(point1, point2);

		if (a.getL() == 1) {
			switch (tau) {
				case 0:
					if (a.geti() == 1) {
						double val1 = x * x * y / (R * R * R * Rxz) +
								x * x * y / (R * Rxz * Rxz * Rxz) -
								y / (R * Rxz);
						double val2 = -x * z / (Rxz * Rxz * Rxz);
						double val3 = x * x / (R * R * R) - 1 / R;
						return new double[]{val1, val2, val3};
					}
					else if (a.getj() == 1) {
						double val1 = x / (R * Rxz) - x * Rxz / (R * R * R);
						double val3 = x * y / (R * R * R);
						return new double[]{val1, 0, val3};
					}
					else if (a.getk() == 1) {
						double val1 = x * y * z / (R * R * R * Rxz) +
								x * y * z / (R * Rxz * Rxz * Rxz);
						double val2 = x * x / (Rxz * Rxz * Rxz) - 1 / Rxz;
						double val3 = x * z / (R * R * R);
						return new double[]{val1, val2, val3};
					}
				case 1:
					if (a.geti() == 1) {
						double val1 =
								x * y * y / (R * R * R * Rxz) - x / (R * Rxz);
						double val3 = x * y / (R * R * R);
						return new double[]{val1, 0, val3};
					}
					else if (a.getj() == 1) {
						double val1 = -y * Rxz / (R * R * R);
						double val3 = y * y / (R * R * R) - 1 / R;
						return new double[]{val1, 0, val3};
					}
					else if (a.getk() == 1) {
						double val1 =
								y * y * z / (R * R * R * Rxz) - z / (R * Rxz);
						double val3 = y * z / (R * R * R);
						return new double[]{val1, 0, val3};
					}
				case 2:
					if (a.geti() == 1) {
						double val1 = x * y * z / (R * R * R * Rxz) +
								x * y * z / (R * Rxz * Rxz * Rxz);
						double val2 = 1 / Rxz - z * z / (Rxz * Rxz * Rxz);
						double val3 = x * z / (R * R * R);
						return new double[]{val1, val2, val3};
					}
					else if (a.getj() == 1) {
						double val1 = z / (R * Rxz) - z * Rxz / (R * R * R);
						double val3 = z * y / (R * R * R);
						return new double[]{val1, 0, val3};
					}
					else if (a.getk() == 1) {
						double val1 = y * z * z / (R * R * R * Rxz) +
								y * z * z / (R * Rxz * Rxz * Rxz) -
								y / (R * Rxz);
						double val2 = x * z / (Rxz * Rxz * Rxz);
						double val3 = z * z / (R * R * R) - 1 / R;
						return new double[]{val1, val2, val3};
					}
			}
		}

		return null;

	}

	public static double gradient(NDDOAtom[] atoms, SolutionR soln,
								  int atomnum,
								  int tau) {

		SimpleMatrix densitymatrix = soln.densityMatrix();


		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.orbsOfAtom;

		int[] atomnumber = soln.atomOfOrb;

		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {


								sum -= atoms[a].getAtomProperties().getQ() *
										GeometryDerivative
												.getGderiv(orbitals[j],
														orbitals[k],
														atoms[a].s(),
														atoms[a].s(), tau);
							}
						}
					}
					else {
						sum -= atoms[atomnum].getAtomProperties().getQ() *
								GeometryDerivative
										.getGderiv(atoms[atomnum].s(),
												atoms[atomnum].s(),
												orbitals[j],
												orbitals[k], tau);
					}
				}
				else {
					if (atomnumber[j] == atomnum) {
						sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) *
								LCGTO.getSDeriv(orbitals[j], orbitals[k], tau);
					}
					else if (atomnumber[k] == atomnum) {
						sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) *
								LCGTO.getSDeriv(orbitals[k], orbitals[j], tau);
					}
				}

				H.set(j, k, sum);
				H.set(k, j, sum);
			}
		}

		SimpleMatrix G = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {

								for (int l : index[a]) {
									if (l > -1) {
										for (int m : index[a]) {
											if (m > -1) {
												sum += densitymatrix.get(l,
														m) *
														GeometryDerivative
																.getGderiv(
																		orbitals[j],
																		orbitals[k],
																		orbitals[l],
																		orbitals[m],
																		tau);
											}
										}
									}
								}
							}
						}
					}
					else {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnum]) {
									if (m > -1) {
										sum += densitymatrix.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[l],
																orbitals[m],
																orbitals[j],
																orbitals[k],
																tau);

									}
								}
							}
						}
					}
				}
				else {
					if (atomnumber[j] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[k]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[j],
																orbitals[l],
																orbitals[k],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[j]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[k],
																orbitals[l],
																orbitals[j],
																orbitals[m],
																tau);
									}
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

		//System.out.println ("Electronic gradient (Analytic): " + e);

		for (int j = 0; j < atoms.length; j++) {
			if (j != atomnum) {
				e += atoms[atomnum].crfDeriv(atoms[j], tau);
			}
		}

//        if (Math.abs(e - grad(atoms, soln, atomnum, tau)) > 1E-5) {
//            System.err.println ("oh well, you can't say you weren't
//            expecting that...");
//            System.exit(0);
//        }


		return e;

	}

	public static double gradientUnrestricted(NDDOAtom[] atoms, SolutionU soln,
											  int atomnum, int tau) {

		SimpleMatrix alphadensity = soln.alphaDensity();

		SimpleMatrix betadensity = soln.betaDensity();

		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.orbsOfAtom;

		int[] atomnumber = soln.atomOfOrb;

		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {


								sum -= atoms[a].getAtomProperties().getQ() *
										GeometryDerivative
												.getGderiv(orbitals[j],
														orbitals[k],
														atoms[a].s(),
														atoms[a].s(), tau);
							}
						}
					}
					else {
						sum -= atoms[atomnum].getAtomProperties().getQ() *
								GeometryDerivative
										.getGderiv(atoms[atomnum].s(),
												atoms[atomnum].s(),
												orbitals[j],
												orbitals[k], tau);
					}
				}
				else {
					if (atomnumber[j] == atomnum) {
						sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) *
								LCGTO.getSDeriv(orbitals[j], orbitals[k], tau);
					}
					else if (atomnumber[k] == atomnum) {
						sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) *
								LCGTO.getSDeriv(orbitals[k], orbitals[j], tau);
					}
				}

				H.set(j, k, sum);
				H.set(k, j, sum);
			}
		}

		SimpleMatrix J = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {

								for (int l : index[a]) {
									if (l > -1) {
										for (int m : index[a]) {
											if (m > -1) {
												sum += (alphadensity.get(l,
														m) +
														betadensity.get(l,
																m)) *
														GeometryDerivative
																.getGderiv(
																		orbitals[j],
																		orbitals[k],
																		orbitals[l],
																		orbitals[m],
																		tau);
											}
										}
									}
								}
							}
						}
					}
					else {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnum]) {
									if (m > -1) {
										sum += (alphadensity.get(l, m) +
												betadensity.get(l, m)) *
												GeometryDerivative
														.getGderiv(orbitals[l],
																orbitals[m],
																orbitals[j],
																orbitals[k],
																tau);
									}
								}
							}
						}
					}
				}

				J.set(j, k, sum);
				J.set(k, j, sum);

			}
		}

		SimpleMatrix Ka = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;

				if (atomnumber[j] != atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[k]]) {
									if (m > -1) {
										sum -= alphadensity.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[j],
																orbitals[l],
																orbitals[k],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[j]]) {
									if (m > -1) {
										sum -= alphadensity.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[k],
																orbitals[l],
																orbitals[j],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
				}

				Ka.set(j, k, sum);
				Ka.set(k, j, sum);

			}
		}

		SimpleMatrix Kb = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;

				if (atomnumber[j] != atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[k]]) {
									if (m > -1) {
										sum -= betadensity.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[j],
																orbitals[l],
																orbitals[k],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[j]]) {
									if (m > -1) {
										sum -= betadensity.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[k],
																orbitals[l],
																orbitals[j],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
				}

				Kb.set(j, k, sum);
				Kb.set(k, j, sum);

			}
		}

		SimpleMatrix Fa = H.plus(J).plus(Ka);
		SimpleMatrix Fb = H.plus(J).plus(Kb);

		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * alphadensity.get(j, k) *
						(H.get(j, k) + Fa.get(j, k));
				e += 0.5 * betadensity.get(j, k) * (H.get(j, k) + Fb.get(j,
						k));
			}
		}


		for (int j = 0; j < atoms.length; j++) {
			if (j != atomnum) {
				e += atoms[atomnum].crfDeriv(atoms[j], tau);
			}
		}

//        if (Math.abs(e - grad(atoms, soln, atomnum, tau)) > 1E-5) {
//            System.err.println ("oh well, you can't say you weren't
//            expecting that...");
//            System.exit(0);
//        }

		return e;
	}

	public static SimpleMatrix[][] gradientRoutine(SolutionR soln) {

		SimpleMatrix[] fockderivatives = new SimpleMatrix[soln.atoms.length * 3];

		SimpleMatrix grad = new SimpleMatrix(soln.atoms.length * 3, 1);

		int count = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			for (int tau = 0; tau < 3; tau++) {
				SimpleMatrix[] matrices = staticderivs(soln.atoms, soln, a, tau);
				fockderivatives[count] = matrices[1];
				double sum = 0;

				for (int i = 0; i < matrices[1].numRows(); i++) {
					for (int j = 0; j < matrices[1].numRows(); j++) {
						sum += 0.5 * soln.densityMatrix().get(i, j) *
								(matrices[0].get(i, j) + matrices[1].get(i,
										j));
					}
				}

				for (int j = 0; j < soln.atoms.length; j++) {
					if (j != a) {
						sum += soln.atoms[a].crfDeriv(soln.atoms[j], tau);
					}
				}

				grad.set(count, 0, sum);

				count++;

			}
		}

		return new SimpleMatrix[][]{new SimpleMatrix[]{grad}, fockderivatives};
	}

	public static SimpleMatrix[][] gradientRoutine(SolutionU soln) {

		SimpleMatrix[] fockderivativesalpha =
				new SimpleMatrix[soln.atoms.length * 3];

		SimpleMatrix[] fockderivativesbeta =
				new SimpleMatrix[soln.atoms.length * 3];

		SimpleMatrix grad = new SimpleMatrix(soln.atoms.length * 3, 1);

		int count = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			for (int tau = 0; tau < 3; tau++) {
				SimpleMatrix[] matrices = staticderivs(soln.atoms, soln, a, tau);
				fockderivativesalpha[count] = matrices[1];
				fockderivativesbeta[count] = matrices[2];
				double sum = 0;

				for (int i = 0; i < matrices[1].numRows(); i++) {
					for (int j = 0; j < matrices[1].numRows(); j++) {
						sum += 0.5 * soln.densityMatrix().get(i, j) *
								(matrices[0].get(i, j) + matrices[1].get(i,
										j));
						sum += 0.5 * soln.densityMatrix().get(i, j) *
								(matrices[0].get(i, j) + matrices[2].get(i,
										j));
					}
				}

				for (int j = 0; j < soln.atoms.length; j++) {
					if (j != a) {
						sum += soln.atoms[a].crfDeriv(soln.atoms[j], tau);
					}
				}

				grad.set(count, 0, sum);

				count++;

			}
		}

		return new SimpleMatrix[][]{new SimpleMatrix[]{grad},
				fockderivativesalpha, fockderivativesbeta};
	}

	public static double grad(SolutionR soln, int atomnum, int tau) {
		double e = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			if (a != atomnum) {
				e += Ederiv(atomnum, a, soln.orbsOfAtom, soln.densityMatrix(),
						soln.atoms, soln.orbitals, tau);
				e += soln.atoms[atomnum].crfDeriv(soln.atoms[a], tau);
			}
		}

		return e;
	}

	public static double grad(SolutionU soln, int atomnum, int tau) {

		double e = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			if (a != atomnum) {
				e += Ederiv(atomnum, a, soln.orbsOfAtom, soln.alphaDensity(),
						soln.betaDensity(), soln.atoms, soln.orbitals, tau);
				e += soln.atoms[atomnum].crfDeriv(soln.atoms[a], tau);
			}
		}

		return e;

	}

	private static double Ederiv(int atomnum1, int atomnum2, int[][] index,
								 SimpleMatrix densityMatrix, NDDOAtom[] atoms,
								 NDDO6G[] orbitals, int tau) {

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += densityMatrix.get(i, j) *
							atoms[atomnum2]
									.Vderiv(orbitals[i], orbitals[j], tau);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e -= densityMatrix.get(k, l) *
							atoms[atomnum1]
									.Vderiv(orbitals[k], orbitals[l], tau);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * densityMatrix.get(i, k) *
							NDDO6G.betaderiv(orbitals[i], orbitals[k], tau);
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
									* GeometryDerivative
									.getGderiv(orbitals[i], orbitals[j],
											orbitals[k], orbitals[l], tau);
						}
					}
				}
			}
		}

		return e;


	}

	private static double Ederiv(int atomnum1, int atomnum2, int[][] index,
								 SimpleMatrix alphaDensity,
								 SimpleMatrix betaDensity, NDDOAtom[] atoms,
								 NDDO6G[] orbitals, int tau) {

		double e = 0;


		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += (alphaDensity.get(i, j) + betaDensity.get(i, j)) *
							atoms[atomnum2]
									.Vderiv(orbitals[i], orbitals[j], tau);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e -= (alphaDensity.get(k, l) + betaDensity.get(k, l)) *
							atoms[atomnum1]
									.Vderiv(orbitals[k], orbitals[l], tau);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * (alphaDensity.get(i, k) + betaDensity.get(i, k)) *
							NDDO6G.betaderiv(orbitals[i], orbitals[k], tau);
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
									* GeometryDerivative
									.getGderiv(orbitals[i], orbitals[j],
											orbitals[k], orbitals[l], tau);
						}
					}
				}
			}
		}

		return e;


	}

	public static SimpleMatrix[] staticderivs(NDDOAtom[] atoms, SolutionR soln,
											  int atomnum, int tau) {

		SimpleMatrix densitymatrix = soln.densityMatrix();


		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.orbsOfAtom;

		int[] atomnumber = soln.atomOfOrb;

		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {


								sum -= atoms[a].getAtomProperties().getQ() *
										GeometryDerivative
												.getGderiv(orbitals[j],
														orbitals[k],
														atoms[a].s(),
														atoms[a].s(), tau);
							}
						}
					}
					else {
						sum -= atoms[atomnum].getAtomProperties().getQ() *
								GeometryDerivative
										.getGderiv(atoms[atomnum].s(),
												atoms[atomnum].s(),
												orbitals[j],
												orbitals[k], tau);
					}
				}
				else {
					if (atomnumber[j] == atomnum) {
						sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) *
								LCGTO.getSDeriv(orbitals[j], orbitals[k], tau);
					}
					else if (atomnumber[k] == atomnum) {
						sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) *
								LCGTO.getSDeriv(orbitals[k], orbitals[j], tau);
					}
				}

				H.set(j, k, sum);
				H.set(k, j, sum);
			}
		}

		SimpleMatrix G = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {

								for (int l : index[a]) {
									if (l > -1) {
										for (int m : index[a]) {
											if (m > -1) {
												sum += densitymatrix.get(l,
														m) *
														GeometryDerivative
																.getGderiv(
																		orbitals[j],
																		orbitals[k],
																		orbitals[l],
																		orbitals[m],
																		tau);
											}
										}
									}
								}
							}
						}
					}
					else {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnum]) {
									if (m > -1) {
										sum += densitymatrix.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[l],
																orbitals[m],
																orbitals[j],
																orbitals[k],
																tau);

									}
								}
							}
						}
					}
				}
				else {
					if (atomnumber[j] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[k]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[j],
																orbitals[l],
																orbitals[k],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[j]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[k],
																orbitals[l],
																orbitals[j],
																orbitals[m],
																tau);
									}
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

		return new SimpleMatrix[]{H, F};

	}

	public static SimpleMatrix[] staticderivs(NDDOAtom[] atoms, SolutionU soln,
											  int atomnum, int tau) {

		SimpleMatrix alphadensity = soln.alphaDensity();

		SimpleMatrix betadensity = soln.betaDensity();

		NDDO6G[] orbitals = soln.orbitals;

		int[][] index = soln.orbsOfAtom;

		int[] atomnumber = soln.atomOfOrb;

		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {


								sum -= atoms[a].getAtomProperties().getQ() *
										GeometryDerivative
												.getGderiv(orbitals[j],
														orbitals[k],
														atoms[a].s(),
														atoms[a].s(), tau);
							}
						}
					}
					else {
						sum -= atoms[atomnum].getAtomProperties().getQ() *
								GeometryDerivative
										.getGderiv(atoms[atomnum].s(),
												atoms[atomnum].s(),
												orbitals[j],
												orbitals[k], tau);
					}
				}
				else {
					if (atomnumber[j] == atomnum) {
						sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) *
								LCGTO.getSDeriv(orbitals[j], orbitals[k], tau);
					}
					else if (atomnumber[k] == atomnum) {
						sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) *
								LCGTO.getSDeriv(orbitals[k], orbitals[j], tau);
					}
				}

				H.set(j, k, sum);
				H.set(k, j, sum);
			}
		}

		SimpleMatrix J = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {

								for (int l : index[a]) {
									if (l > -1) {
										for (int m : index[a]) {
											if (m > -1) {
												sum += (alphadensity.get(l,
														m) +
														betadensity.get(l,
																m)) *
														GeometryDerivative
																.getGderiv(
																		orbitals[j],
																		orbitals[k],
																		orbitals[l],
																		orbitals[m],
																		tau);
											}
										}
									}
								}
							}
						}
					}
					else {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnum]) {
									if (m > -1) {
										sum += (alphadensity.get(l, m) +
												betadensity.get(l, m)) *
												GeometryDerivative
														.getGderiv(orbitals[l],
																orbitals[m],
																orbitals[j],
																orbitals[k],
																tau);
									}
								}
							}
						}
					}
				}

				J.set(j, k, sum);
				J.set(k, j, sum);

			}
		}

		SimpleMatrix Ka = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;

				if (atomnumber[j] != atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[k]]) {
									if (m > -1) {
										sum -= alphadensity.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[j],
																orbitals[l],
																orbitals[k],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[j]]) {
									if (m > -1) {
										sum -= alphadensity.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[k],
																orbitals[l],
																orbitals[j],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
				}

				Ka.set(j, k, sum);
				Ka.set(k, j, sum);

			}
		}

		SimpleMatrix Kb = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;

				if (atomnumber[j] != atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[k]]) {
									if (m > -1) {
										sum -= betadensity.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[j],
																orbitals[l],
																orbitals[k],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l : index[atomnum]) {
							if (l > -1) {
								for (int m : index[atomnumber[j]]) {
									if (m > -1) {
										sum -= betadensity.get(l, m) *
												GeometryDerivative
														.getGderiv(orbitals[k],
																orbitals[l],
																orbitals[j],
																orbitals[m],
																tau);
									}
								}
							}
						}
					}
				}

				Kb.set(j, k, sum);
				Kb.set(k, j, sum);

			}
		}

		SimpleMatrix Fa = H.plus(J).plus(Ka);
		SimpleMatrix Fb = H.plus(J).plus(Kb);

		return new SimpleMatrix[]{H, Fa, Fb};

	}

	public static SimpleMatrix densitymatrixderivfinite(NDDOAtom[] atoms,
														SolutionR soln,
														int atomnum, int tau) {
		SimpleMatrix orig = soln.densityMatrix();

		NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum, tau);

		SimpleMatrix perturbed =
				soln.withNewAtoms(newatoms).densityMatrix();

		return perturbed.minus(orig).scale(1E7);
	}

	public static SimpleMatrix[] densitymatrixderivfinite(NDDOAtom[] atoms,
														  SolutionU soln,
														  int atomnum,
														  int tau) {
		SimpleMatrix aorig = soln.alphaDensity();
		SimpleMatrix borig = soln.betaDensity();

		NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum, tau);
		SolutionU newsoln = (SolutionU) soln.withNewAtoms(newatoms);

		SimpleMatrix aperturbed = newsoln.alphaDensity();
		SimpleMatrix bperturbed = newsoln.betaDensity();

		return new SimpleMatrix[]{aperturbed.minus(aorig).scale(1E7),
				bperturbed.minus(borig).scale(1E7)};
	}
}
