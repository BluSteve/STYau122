package nddo.param;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.NDDOParams;
import nddo.solution.SolutionR;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static nddo.State.nom;
import static nddo.geometry.GeometrySecondDerivative.computeResponseVectorsPople;

public class ParamDerivative {
	@Deprecated
	public static double crfderivfinite(NDDOAtom A, NDDOAtom B, int num) {

		double initial = A.crf(B);

		if (num == 0 || num == 2) {

			try {
				NDDOParams params = A.getParams();
				params.modifyParam(0, Constants.LAMBDA);

				A = A.withNewParams(params);
			} catch (Exception e) {
				e.printStackTrace();
//				System.exit(0);
			}
		}
		if (num == 1 || num == 2) {
			try {
				NDDOParams params = B.getParams();
				params.modifyParam(0, Constants.LAMBDA);

				B = B.withNewParams(params);
			} catch (Exception e) {
				e.printStackTrace();
//				System.exit(0);
			}
		}

		double finalval = A.crf(B);

		return (finalval - initial) / Constants.LAMBDA;


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


		NDDOOrbital[] orbitals = soln.orbitals;

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
					double Huk = nom.betaparamderiv(orbitals[j],
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
												nom.getGderiv(orbitals[j],
														orbitals[k],
														orbitals[l],
														orbitals[m],
														getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]],
																Z),
														type);

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
													nom.getGderiv(
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


		NDDOOrbital[] orbitals = soln.orbitals;

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
					double Huk = nom.betaparamderiv(orbitals[j],
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
												nom.getGderiv(orbitals[j],
														orbitals[k],
														orbitals[l],
														orbitals[m],
														getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]],
																Z),
														type);

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
													nom
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

		return H.copy().plus(G);
	}

	private static SimpleMatrix zetaHderivstatic(NDDOAtom[] atoms,
												 SolutionR soln, int Z,
												 int type) {


		NDDOOrbital[] orbitals = soln.orbitals;

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
					double Huk = nom.betaparamderiv(orbitals[j],
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


		NDDOOrbital[] orbitals = soln.orbitals;

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
												nom
														.getGderiv(orbitals[j],
																orbitals[k],
																orbitals[l],
																orbitals[m],
																getNum(atomicnumbers[atomNumber[j]],
																		atomicnumbers[atomNumber[l]],
																		Z),
																type);

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
													nom
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


		NDDOOrbital[] orbitals = soln.orbitals;

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


		NDDOOrbital[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;


		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] != atomNumber[k]) {
					double sum = betaderivsum(atomicnumbers[atomNumber[j]],
							atomicnumbers[atomNumber[k]], Z,
							orbitals[j].getL(),
							orbitals[k].getL(), type);

					if (sum != 0) {
						double H = nom.betaparambetaderiv(orbitals[j], orbitals[k], sum);
						e += 2 * densitymatrix.get(j, k) * H;
					}
				}
			}
		}

		return e / 4.3363E-2;

	}

	private static SimpleMatrix betafockderivstatic(SolutionR soln, int Z,
													int type) {

		SimpleMatrix F =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		NDDOOrbital[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] != atomNumber[k]) {
					double sum = betaderivsum(atomicnumbers[atomNumber[j]],
							atomicnumbers[atomNumber[k]], Z,
							orbitals[j].getL(),
							orbitals[k].getL(), type);

					if (sum != 0) {
						double H = nom.betaparambetaderiv(orbitals[j], orbitals[k], sum);
						F.set(j, k, H);
						F.set(k, j, H);
					}
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
				sum += soln.atoms[i].crfAlphaDeriv(soln.atoms[j],
						getNum(soln.atomicNumbers[i], soln.atomicNumbers[j],
								Z));
			}
		}

		return sum / 4.3363E-2;
	}

	private static double betaderivsum(int Z1, int Z2, int Z, int L1, int L2,
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
		double[] Darr = new double[nonv];
		double[] Dinvarr = new double[nonv];

		int counter = 0;
		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = -soln.E.get(i) + soln.E.get(NOcc + j);

				Darr[counter] = Math.pow(e, -0.5);
				Dinvarr[counter] = Math.pow(e, 0.5);

				counter++;
			}
		}

		// convert AO to MO basis
		SimpleMatrix F = new SimpleMatrix(nonv, length);
		for (int a = 0; a < length; a++) {
			SimpleMatrix f = new SimpleMatrix(nonv, 1);

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

					element /= soln.E.get(j + NOcc) - soln.E.get(i);

					f.set(count, 0, element);

					count++;
				}
			}

			CommonOps_DDRM.multRows(Darr, f.getDDRM());
			barray[a] = f;
			Farray[a] = barray[a].copy();
			F.setColumn(a, 0, barray[a].getDDRM().data);
		}

		// main loop
		int[] iterable = new int[length];

		// 0: B, 1: Bt, 2: Bn, 3: P, 4: BmP
		List<SimpleMatrix[]> prevs = new ArrayList<>();
		List<Double> dots = new ArrayList<>();

		double[] oldrMags = new double[rarray.length];
		Arrays.fill(oldrMags, 1);

		bigLoop:
		while (Utils.numIterable(iterable) > 0) {
			// orthogonalize barray
			for (int i = 0; i < barray.length; i++) {
				for (int j = 0; j < i; j++) {
					barray[i].plusi(barray[i].dot(barray[j]) /
							barray[j].dot(barray[j]), barray[j].negativei());
				}
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix[] prev = new SimpleMatrix[5];
				prev[0] = barray[i]; // original barray object here
				prev[1] = barray[i].transpose();
				prev[2] = barray[i].negative();
				dots.add(barray[i].dot(barray[i]));

				// parray[i] stays the same object throughout
				SimpleMatrix bc = barray[i].copy();
				CommonOps_DDRM.multRows(Dinvarr, bc.getDDRM());
				SimpleMatrix crv = computeResponseVectorsPople(bc, soln);
				CommonOps_DDRM.multRows(Darr, crv.getDDRM());
				parray[i] = crv;

				prev[3] = parray[i];
				prev[4] = barray[i].minus(parray[i]);

				prevs.add(prev);
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix newb = parray[i].copy();

				// orthogonalize against all previous Bs
				for (int j = 0; j < prevs.size(); j++) {
					SimpleMatrix[] prev = prevs.get(j);
					SimpleMatrix transpose = prev[1];
					double num = transpose.mult(parray[i]).get(0) /
							dots.get(j);

					newb.plusi(num, prev[2]);
				}

				barray[i] = newb; // new barray object created
			}

			SimpleMatrix Bt = new SimpleMatrix(prevs.size(), nonv);
			SimpleMatrix P = new SimpleMatrix(nonv, prevs.size());
			for (int i = 0; i < prevs.size(); i++) {
				Bt.setRow(i, 0, prevs.get(i)[0].getDDRM().data);
				P.setColumn(i, 0, prevs.get(i)[4].getDDRM().data);
			}

			SimpleMatrix rhs = Bt.mult(F);
			SimpleMatrix lhs = Bt.mult(P);

			// alpha dimensions are prevBs x length
			SimpleMatrix alpha;
			try {
				alpha = lhs.solve(rhs);
			} catch (SingularMatrixException e) {
				alpha = SimpleMatrix.ones(lhs.numCols(), rhs.numCols());
			}

			// reset r and x array
			for (int a = 0; a < length; a++) {
				rarray[a] = new SimpleMatrix(nonv, 1);
				xarray[a] = new SimpleMatrix(nonv, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {
					// B with tilde
					rarray[j].plusi(alpha.get(i, j), prevs.get(i)[4]);
					xarray[j].plusi(alpha.get(i, j), prevs.get(i)[0]);
				}
			}

			int xarrayHoldNN = Utils.numNotNull(xarrayHold);
			for (int j = 0; j < alpha.numCols(); j++) {
				rarray[j].minusi(Farray[j]);
				CommonOps_DDRM.multRows(Dinvarr, xarray[j].getDDRM());

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
				alpha = SimpleMatrix.ones(solver.numCols(), rhsvec.numCols());
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
				beta = SimpleMatrix.ones(solver.numCols(), rhsvec.numCols());
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
			if (soln.E.get(NOcc - 1) != soln.E.get(j)) {
				x.set(count1, 0,
						-element / (soln.E.get(NOcc - 1) - soln.E.get(j)));
			}
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
					D1deriv = atoms[i].D1Deriv(paramnum - 5);
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

			if (index[j].length > 1) { // exclude hydrogen
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

	private static double D1Derivfinite(NDDOAtom a, int type) {

		double D1 = a.D1;


		NDDOParams params = a.getParams();
		params.modifyParam(5 + type, Constants.LAMBDA);

		NDDOAtom a2 = a.withNewParams(params);

		double D1perturbed = a2.D1;

		return (D1perturbed - D1) / Constants.LAMBDA;

	}

	private static double D2Derivfinite(NDDOAtom a, int type) {

		double D2 = a.D2;


		NDDOParams params = a.getParams();
		params.modifyParam(5 + type, Constants.LAMBDA);

		NDDOAtom a2 = a.withNewParams(params);

		double D2perturbed = a2.D2;

		return (D2perturbed - D2) / Constants.LAMBDA;

	}

	private static double p1Derivfinite(NDDOAtom a, int type) {

		double p1 = a.p1;

		NDDOParams params = a.getParams();
		params.modifyParam(5 + type, Constants.LAMBDA);

		NDDOAtom a2 = a.withNewParams(params);

		double p1perturbed = a2.p1;

		return (p1perturbed - p1) / Constants.LAMBDA;

	}

	private static double p2Derivfinite(NDDOAtom a, int type) {

		double p2 = a.p2;

		NDDOParams params = a.getParams();
		params.modifyParam(5 + type, Constants.LAMBDA);

		NDDOAtom a2 = a.withNewParams(params);

		double p2perturbed = a2.p2;

		return (p2perturbed - p2) / Constants.LAMBDA;

	}
}
