package nddo.param;

import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.geometry.GeometrySecondDerivative;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.apache.commons.lang3.time.StopWatch;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import tools.Pow;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static nddo.State.nom;
import static nddo.math.PopleThiel.computeResponseVectorsPople;
import static tools.Utils.mag;

public class ParamDerivative {
	public static double HFDeriv(Solution soln, int Z, int paramnum) {

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
			if (soln instanceof SolutionR) {
				return zetaHfderiv((SolutionR) soln, Z, paramnum - 5);
			}
			else if (soln instanceof SolutionU) {
				return zetaHfderiv((SolutionU) soln, Z, paramnum - 5);
			}
		}
		if (paramnum == 7) {
			return eisolHfderiv(soln.atoms, Z);
		}

		System.err.println("oh no! This isn't MNDO!");
		return 0;

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

	public static SimpleMatrix[][] MNDOStaticMatrixDeriv(SolutionU soln, int Z, int firstParamIndex) {
		NDDOAtom[] atoms = soln.atoms;
		SimpleMatrix[] HDerivs = new SimpleMatrix[8];
		SimpleMatrix[] FaDerivs = new SimpleMatrix[8];
		SimpleMatrix[] FbDerivs = new SimpleMatrix[8];


		if (firstParamIndex <= 1) HDerivs[1] = betafockderivstatic(soln, Z, 0);
		if (firstParamIndex <= 1) FaDerivs[1] = HDerivs[1].copy();
		if (firstParamIndex <= 1) FbDerivs[1] = HDerivs[1].copy();
		if (firstParamIndex <= 3) HDerivs[3] = uxxfockderivstatic(soln, Z, 0);
		if (firstParamIndex <= 3) FaDerivs[3] = HDerivs[3].copy();
		if (firstParamIndex <= 3) FbDerivs[3] = HDerivs[3].copy();

		if (firstParamIndex <= 5)
			HDerivs[5] = zetaHderivstatic(atoms, soln, Z, 0);
		if (firstParamIndex <= 5) {

			SimpleMatrix[] matrices = zetaGderivstatic(atoms, soln, Z, 0);
			FaDerivs[5] = HDerivs[5].copy().plus(matrices[0]);
			FbDerivs[5] = HDerivs[5].copy().plus(matrices[1]);

		}

		if (Z != 1) {
			if (firstParamIndex <= 2)
				HDerivs[2] = betafockderivstatic(soln, Z, 1);
			if (firstParamIndex <= 2) FaDerivs[2] = HDerivs[2].copy();
			if (firstParamIndex <= 2) FbDerivs[2] = HDerivs[2].copy();
			if (firstParamIndex <= 4)
				HDerivs[4] = uxxfockderivstatic(soln, Z, 1);
			if (firstParamIndex <= 4) FaDerivs[4] = HDerivs[4].copy();
			if (firstParamIndex <= 4) FbDerivs[4] = HDerivs[4].copy();
			if (firstParamIndex <= 6)
				HDerivs[6] = zetaHderivstatic(atoms, soln, Z, 1);
			if (firstParamIndex <= 6) {
				SimpleMatrix[] matrices = zetaGderivstatic(atoms, soln, Z, 1);
				FaDerivs[6] = HDerivs[6].copy().plus(matrices[0]);
				FbDerivs[6] = HDerivs[6].copy().plus(matrices[1]);

			}
		}
		return new SimpleMatrix[][]{HDerivs, FaDerivs, FbDerivs};
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
					for (int an = 0; an < soln.atoms.length; an++)
						if (atomNumber[j] != an) Huv += soln.atoms[an].Vpd(orbitals[j], orbitals[k],
								getNum(atomicnumbers[an], atomicnumbers[atomNumber[j]], Z), type);
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { /* case 3*/
					double Huk = nom.Hzetapd(orbitals[j], orbitals[k],
							getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]], Z), type);
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
					for (int l : missingIndex[atomNumber[j]])
						if (l > -1) for (int m : missingIndex[atomNumber[j]])
							if (m > -1 && atomNumber[l] == atomNumber[m]) sum += soln.densityMatrix().get(l, m) *
									nom.Gpd(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
											getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[l]], Z),
											type);
				}
				else for (int l : index[atomNumber[j]])
					if (l > -1) for (int m : index[atomNumber[k]])
						if (m > -1) sum += soln.densityMatrix().get(l, m) * (-0.5 *
								nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
										getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]], Z), type));
				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		}

		SimpleMatrix F = H.copy().plus(G);
		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) e += 0.5 * densitymatrix.get(j, k) * (H.get(j, k) + F.get(j, k));
		}
		return e / 4.3363E-2;
	}

	private static double zetaHfderiv(SolutionU soln, int Z, int type) {

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
									.Vpd(orbitals[j], orbitals[k],
											getNum(atomicnumbers[an],
													atomicnumbers[atomNumber[j]],
													Z), type);
						}
					}
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { // case 3
					double Huk = nom.Hzetapd(orbitals[j], orbitals[k],
							getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]], Z), type);
					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}

		SimpleMatrix Ga = new SimpleMatrix(orbitals.length, orbitals.length);

		SimpleMatrix Gb = new SimpleMatrix(orbitals.length, orbitals.length);


		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double suma = 0;
				double sumb = 0;

				if (atomNumber[j] == atomNumber[k]) {
					for (int l : missingIndex[atomNumber[j]]) {
						if (l > -1) {
							for (int m : missingIndex[atomNumber[j]]) {
								if (m > -1) {
									if (atomNumber[l] == atomNumber[m]) {
										suma += soln.densityMatrix().get(l, m) *
												nom.Gpd(orbitals[j], orbitals[k], orbitals[l],
														orbitals[m], getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z), type);
										sumb += soln.densityMatrix().get(l, m) *
												nom.Gpd(orbitals[j], orbitals[k], orbitals[l],
														orbitals[m], getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z), type);

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
									suma -= soln.alphaDensity().get(l, m) *
											nom.Gpd(orbitals[j], orbitals[l], orbitals[k],
													orbitals[m],
													getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]],
															Z), type);

									sumb -= soln.betaDensity().get(l, m) *
											nom.Gpd(orbitals[j], orbitals[l], orbitals[k],
													orbitals[m],
													getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]],
															Z), type);

								}
							}
						}
					}
				}

				Ga.set(j, k, suma);
				Ga.set(k, j, suma);

				Gb.set(j, k, sumb);
				Gb.set(k, j, sumb);
			}
		}

		SimpleMatrix Fa = H.copy().plus(Ga);
		SimpleMatrix Fb = H.copy().plus(Gb);


		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * soln.alphaDensity().get(j, k) * (H.get(j, k) + Fa.get(j, k));
				e += 0.5 * soln.betaDensity().get(j, k) * (H.get(j, k) + Fb.get(j, k));

			}
		}

		return e / 4.3363E-2;

	}

	private static SimpleMatrix zetaHderivstatic(NDDOAtom[] atoms,
												 Solution soln, int Z,
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
									.Vpd(orbitals[j], orbitals[k],
											getNum(atomicnumbers[an],
													atomicnumbers[atomNumber[j]],
													Z), type);
						}
					}
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { // case 3
					double Huk = nom.Hzetapd(orbitals[j],
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
														.Gpd(orbitals[j],
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
															.Gpd(
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

	private static SimpleMatrix[] zetaGderivstatic(NDDOAtom[] atoms, SolutionU soln, int Z, int type) {
		NDDOOrbital[] orbitals = soln.orbitals;

		int[][] index = soln.orbsOfAtom;

		int[][] missingIndex = soln.missingOfAtom;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		SimpleMatrix Ga = new SimpleMatrix(orbitals.length, orbitals.length);

		SimpleMatrix Gb = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double suma = 0;
				double sumb = 0;

				if (atomNumber[j] == atomNumber[k]) {
					for (int l : missingIndex[atomNumber[j]]) {
						if (l > -1) {
							for (int m : missingIndex[atomNumber[j]]) {
								if (m > -1) {
									if (atomNumber[l] == atomNumber[m]) {
										suma += soln.densityMatrix().get(l, m) *
												State.nom.Gpd(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
														getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z), type);
										sumb += soln.densityMatrix().get(l, m) *
												State.nom.Gpd(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
														getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z), type);

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
									suma -= soln.alphaDensity().get(l, m) *
											State.nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
													getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]],
															Z), type);

									sumb -= soln.betaDensity().get(l, m) *
											State.nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
													getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]],
															Z), type);

								}
							}
						}
					}
				}

				Ga.set(j, k, suma);
				Ga.set(k, j, suma);

				Gb.set(j, k, sumb);
				Gb.set(k, j, sumb);
			}
		}

		return new SimpleMatrix[]{Ga, Gb};
	}


	private static double uxxHfderiv(Solution soln, int Z, int type) {

		SimpleMatrix densitymatrix = soln.densityMatrix();


		NDDOOrbital[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {

			if (atomicnumbers[atomNumber[j]] == Z && orbitals[j].getL() == type) {
				e += densitymatrix.get(j, j);
			}
		}

		return e / 4.3363E-2;

	}

	private static SimpleMatrix uxxfockderivstatic(Solution soln, int Z,
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

	private static double betaHfderiv(Solution soln, int Z, int type) {

		SimpleMatrix densitymatrix = soln.densityMatrix();


		NDDOOrbital[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;


		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] != atomNumber[k]) {
					double H = nom.Hbetapd(orbitals[j], orbitals[k],
							getNumBeta(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]],
									Z, orbitals[j].getL(), orbitals[k].getL(), type));
					e += 2 * densitymatrix.get(j, k) * H;
				}
			}
		}

		return e / 4.3363E-2;
	}

	private static SimpleMatrix betafockderivstatic(Solution soln, int Z,
													int type) {

		SimpleMatrix F =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		NDDOOrbital[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] != atomNumber[k]) {
					double H = nom.Hbetapd(orbitals[j], orbitals[k],
							getNumBeta(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]],
									Z, orbitals[j].getL(), orbitals[k].getL(), type));
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

	private static double alphaHfderiv(Solution soln, int Z) {

		double sum = 0;

		for (int i = 0; i < soln.atoms.length; i++) {
			for (int j = i + 1; j < soln.atoms.length; j++) {
				sum += soln.atoms[i].crfalphapd(soln.atoms[j],
						getNum(soln.atomicNumbers[i], soln.atomicNumbers[j],
								Z));
			}
		}

		return sum / 4.3363E-2;
	}

	static int getNum(int Z1, int Z2, int Z) {
		int num = 0;

		if (Z1 == Z) {
			num += 1;
		}

		if (Z2 == Z) {
			num += 2;
		}

		return num - 1;
	}

	static int getNumBeta(int Z1, int Z2, int Z, int L1, int L2, int L) {
		int num = 0;

		if (Z1 == Z && L1 == L) num += 1;
		if (Z2 == Z && L2 == L) num += 2;

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

	public static SimpleMatrix[] responseMatrices(SolutionU soln, SimpleMatrix[] densityderivs) {

		SimpleMatrix densityderivalpha = densityderivs[0];

		SimpleMatrix densityderivbeta = densityderivs[1];


		SimpleMatrix Jderiv = new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		SimpleMatrix Kaderiv = new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		SimpleMatrix Kbderiv = new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);


		int Jcount = 0;
		int Kcount = 0;

		//construct J matrix

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = j; k < soln.orbitals.length; k++) {
				double val = 0;
				if (j == k) {

					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							val += (densityderivalpha.get(l, l) +
									densityderivbeta.get(l, l)) *
									soln.integralArrayCoulomb[Jcount];
							Jcount++;
						}
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
										val += (densityderivalpha.get(l, m) +
												densityderivbeta.get(l, m)) *
												soln.integralArrayCoulomb[Jcount];
										Jcount++;
									}
								}

							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += (densityderivalpha.get(j, k) +
							densityderivbeta.get(j, k)) *
							soln.integralArrayCoulomb[Jcount];
					Jcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
										val += (densityderivalpha.get(l, m) +
												densityderivbeta.get(l, m)) *
												soln.integralArrayCoulomb[Jcount];
										Jcount++;
									}
								}

							}
						}
					}
				}


				Jderiv.set(j, k, val);
				Jderiv.set(k, j, val);
			}
		}

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = j; k < soln.orbitals.length; k++) {
				double vala = 0;
				double valb = 0;
				if (j == k) {

					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							vala += densityderivalpha.get(l, l) *
									soln.integralArrayExchange[Kcount];
							valb += densityderivbeta.get(l, l) *
									soln.integralArrayExchange[Kcount];
							Kcount++;
						}
					}

				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					vala += densityderivalpha.get(j, k) *
							soln.integralArrayExchange[Kcount];
					valb += densityderivbeta.get(j, k) *
							soln.integralArrayExchange[Kcount];
					Kcount++;

				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : soln.orbsOfAtom[soln.atomOfOrb[k]]) {
								if (m > -1) {
									vala += densityderivalpha.get(l, m) *
											soln.integralArrayExchange[Kcount];
									valb += densityderivbeta.get(l, m) *
											soln.integralArrayExchange[Kcount];
									Kcount++;
								}
							}
						}
					}
				}

				Kaderiv.set(j, k, vala);
				Kaderiv.set(k, j, vala);
				Kbderiv.set(j, k, valb);
				Kbderiv.set(k, j, valb);
			}
		}

		SimpleMatrix responsealpha = Jderiv.plus(Kaderiv);

		SimpleMatrix responsebeta = Jderiv.plus(Kbderiv);


		return new SimpleMatrix[]{responsealpha, responsebeta};
	}

	private static SimpleMatrix computeResponseVectorsLimited(SimpleMatrix x,
															  SolutionR soln) {//todo
		// copylicate from GeometrySecondDerivative

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
						sum -= 2 * (soln.Ct.get(i, u) * soln.Ct.get(j + NOcc,
								v) +
								soln.Ct.get(j + NOcc, u) * soln.Ct.get(i, v)) *
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
						element += soln.Ct.get(i, u) * soln.Ct.get(j + NOcc, v) *
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

	public static SimpleMatrix[] xArrayPople(SolutionR soln,
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

				Darr[counter] = Pow.pow(e, -0.5);
				Dinvarr[counter] = Pow.pow(e, 0.5);

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
							element += soln.Ct.get(i, u) *
									soln.Ct.get(j + NOcc, v) *
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
				SimpleMatrix crv = computeResponseVectorsPople(soln, bc);
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
							return xArrayThiel(soln, fockderivstatic);
						}
						if (mag != mag) {
							soln.getRm().getLogger()
									.error("Pople algorithm fails; " +
											"reverting to Thiel algorithm...");
							return xArrayThiel(soln, fockderivstatic);
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

	private static SimpleMatrix[] xArrayThiel(SolutionR soln,
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

				arrpreconditioner[counter] = Pow.pow(e, -0.5);

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
							element += soln.Ct.get(i, u) * soln.Ct.get(j + NOcc,
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
					p.add(D.mult(computeResponseVectorsLimited(dirs[i], soln)));
				}
			}

			SimpleMatrix solver = new SimpleMatrix(p.size(), p.size());
			SimpleMatrix rhsvec = new SimpleMatrix(p.size(), rarray.length);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					double[] arrrhs = new double[p.size()];

					for (int i = 0; i < arrrhs.length; i++) {
						arrrhs[i] = 2 * rarray[a].transpose().mult(d.get(i)).get(0, 0);

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

	public static SimpleMatrix[] xArrayThiel(SolutionU soln, SimpleMatrix[] fockderivstaticalpha,
											 SimpleMatrix[] fockderivstaticbeta) {

		StopWatch sw = new StopWatch();
		sw.start();

		int NOccAlpha = soln.getRm().nOccAlpha;
		int NOccBeta = soln.getRm().nOccBeta;

		int NVirtAlpha = soln.getRm().nVirtAlpha;
		int NVirtBeta = soln.getRm().nVirtBeta;

		SimpleMatrix[] xarray = new SimpleMatrix[fockderivstaticalpha.length];

		SimpleMatrix[] rarray = new SimpleMatrix[fockderivstaticalpha.length];

		SimpleMatrix[] dirs = new SimpleMatrix[fockderivstaticalpha.length];

		SimpleMatrix preconditioner = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);

		int counter = 0;

		for (int i = 0; i < NOccAlpha; i++) {
			for (int j = 0; j < NVirtAlpha; j++) {
				double e = -soln.Ea.get(i) + soln.Ea.get(NOccAlpha + j);
				preconditioner.set(counter, Pow.pow(e, -0.5));
				counter++;
			}
		}

		for (int i = 0; i < NOccBeta; i++) {
			for (int j = 0; j < NVirtBeta; j++) {
				double e = -soln.Eb.get(i) + soln.Eb.get(NOccBeta + j);
				preconditioner.set(counter, Pow.pow(e, -0.5));
				counter++;
			}
		}

		SimpleMatrix D = SimpleMatrix.diag(preconditioner.getDDRM().data);

		//SimpleMatrix D = SimpleMatrix.eye(NOcc * NVirt);

		for (int a = 0; a < xarray.length; a++) {
			SimpleMatrix F = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);

			int count1 = 0;

			for (int i = 0; i < NOccAlpha; i++) { // kappa
				for (int j = 0; j < NVirtAlpha; j++) { // i

					double element = 0;

					for (int u = 0; u < soln.orbitals.length; u++) {
						for (int v = 0; v < soln.orbitals.length; v++) {
							element += soln.ca.get(i, u) * soln.ca.get(j + NOccAlpha, v) * fockderivstaticalpha[a].get(u, v);
						}
					}


					F.set(count1, 0, element);

					count1++;
				}
			}

			for (int i = 0; i < NOccBeta; i++) { // kappa
				for (int j = 0; j < NVirtBeta; j++) { // i

					double element = 0;

					for (int u = 0; u < soln.orbitals.length; u++) {
						for (int v = 0; v < soln.orbitals.length; v++) {
							element += soln.cb.get(i, u) * soln.cb.get(j + NOccBeta, v) * fockderivstaticbeta[a].get(u, v);
						}
					}

					F.set(count1, 0, element);

					count1++;
				}
			}

			F = D.mult(F);

			xarray[a] = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);

			rarray[a] = F.copy();

			dirs[a] = F.copy();
		}


		int numit = 0;


		while (Utils.numNotNull(rarray) > 0) {

			numit++;

			ArrayList<SimpleMatrix> d = new ArrayList<>();

			ArrayList<SimpleMatrix> p = new ArrayList<>();

			//System.err.println("It's still running, don't worry: " + Utils.numNotNull(rarray));

			for (int i = 0; i < rarray.length; i++) {

				if (rarray[i] != null) {

					d.add(dirs[i].copy());
					p.add(D.mult(GeometrySecondDerivative.computeResponseVectorsThiel(dirs[i], soln)));
				}


			}


			SimpleMatrix solver = new SimpleMatrix(p.size(), p.size());


			SimpleMatrix rhsvec = new SimpleMatrix(p.size(), rarray.length);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					SimpleMatrix rhs = new SimpleMatrix(p.size(), 1);

					for (int i = 0; i < rhs.numRows(); i++) {
						rhs.set(i, 0, 2 * rarray[a].transpose().mult(d.get(i)).get(0, 0));

					}

					rhsvec.setColumn(a, 0, rhs.getDDRM().data);
				}
			}

			for (int i = 0; i < solver.numRows(); i++) {
				for (int j = i; j < solver.numRows(); j++) {

					double val = p.get(j).transpose().mult(d.get(i)).get(0,
							0) + p.get(i).transpose().mult(d.get(j)).get(0, 0);
					solver.set(i, j, val);
					solver.set(j, i, val);
				}
			}

			SimpleMatrix alpha = solver.solve(rhsvec);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {


					for (int i = 0; i < alpha.numRows(); i++) {
						xarray[a] =
								xarray[a].plus(d.get(i).scale(alpha.get(i, a)));
						rarray[a] =
								rarray[a].minus(p.get(i).scale(alpha.get(i, a)));

					}

					if (mag(rarray[a]) < 1E-10) {//todo change this
						rarray[a] = null;
					}
					else {
						//System.out.println("convergence test: " + mag(rarray[a]));

					}

				}
			}


			solver = new SimpleMatrix(solver.numRows(), solver.numRows());

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					SimpleMatrix rhs = new SimpleMatrix(solver.numRows(), 1);

					for (int i = 0; i < rhs.numRows(); i++) {
						rhs.set(i, 0, -rarray[a].transpose().mult(p.get(i))
								.get(0, 0));

					}

					rhsvec.setColumn(a, 0, rhs.getDDRM().data);
				}
			}


			for (int i = 0; i < solver.numRows(); i++) {
				for (int j = 0; j < solver.numRows(); j++) {
					solver.set(i, j, d.get(j).transpose().mult(p.get(i)).get(0, 0));
				}
			}

			SimpleMatrix beta = solver.solve(rhsvec);
			for (int a = 0; a < rhsvec.numCols(); a++) {

				if (rarray[a] != null) {


					dirs[a] = rarray[a].copy();

					for (int i = 0; i < beta.numRows(); i++) {
						dirs[a] = dirs[a].plus(d.get(i).scale(beta.get(i, a)));
					}
				}
			}


		}
		return xarray;
	}

	public static SimpleMatrix[] densityDerivatives(SolutionU soln, SimpleMatrix xvector) {
		int NOccAlpha = soln.getRm().nOccAlpha;
		int NOccBeta = soln.getRm().nOccBeta;

		int NVirtAlpha = soln.getRm().nVirtAlpha;
		int NVirtBeta = soln.getRm().nVirtBeta;

		SimpleMatrix densityderivalpha = new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		SimpleMatrix densityderivbeta = new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);


		for (int u = 0; u < densityderivalpha.numRows(); u++) {
			for (int v = u; v < densityderivalpha.numCols(); v++) {
				double sum = 0;
				int count = 0;
				for (int i = 0; i < NOccAlpha; i++) {
					for (int j = 0; j < NVirtAlpha; j++) {
						sum -= (soln.ca.get(i, u) * soln.ca.get(j + NOccAlpha, v) +
								soln.ca.get(j + NOccAlpha, u) * soln.ca.get(i, v)) *
								xvector.get(count, 0);
						count++;
					}
				}

				densityderivalpha.set(u, v, sum);
				densityderivalpha.set(v, u, sum);
			}
		}

		for (int u = 0; u < densityderivalpha.numRows(); u++) {
			for (int v = u; v < densityderivalpha.numCols(); v++) {
				double sum = 0;
				int count = NOccAlpha * NVirtAlpha;
				for (int i = 0; i < NOccBeta; i++) {
					for (int j = 0; j < NVirtBeta; j++) {
						sum -= (soln.cb.get(i, u) * soln.cb.get(j + NOccBeta, v) +
								soln.cb.get(j + NOccBeta, u) * soln.cb.get(i, v)) *
								xvector.get(count, 0);
						count++;
					}
				}

				densityderivbeta.set(u, v, sum);
				densityderivbeta.set(v, u, sum);
			}
		}

		return new SimpleMatrix[]{densityderivalpha, densityderivbeta};
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
					element += soln.Ct.get(NOcc - 1, u) * soln.Ct.get(j, v) *
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
						sum -= 2 * (soln.Ct.get(i, u) * soln.Ct.get(j + NOcc,
								v) +
								soln.Ct.get(j + NOcc, u) *
										soln.Ct.get(i, v)) *
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
					sum -= soln.Ct.get(k, u) * x.get(k, 0);
				}
				else if (k >= NOcc) {
					if (k > 0) {
						sum -= soln.Ct.get(k, u) * x.get(k - 1, 0);
					}
				}
			}

			CDeriv.set(0, u, sum);
		}

		return CDeriv;
	}

	// todo rewrite
	public static double MNDOHomoDerivtemp(SolutionR soln, SimpleMatrix coeffDeriv,
										   SimpleMatrix Fderiv) {

		int index = (int) (soln.nElectrons / 2.0) - 1;

		if (index < 0) {
			return 0;
		}

		SimpleMatrix coeff = soln.Ct.extractVector(true, index);

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

	public static double MNDOHomoDerivNew(SolutionR soln, SimpleMatrix x, SimpleMatrix Fderiv) {

		SimpleMatrix Caderiv = soln.C.mult(x);
		SimpleMatrix Eaderiv = Caderiv.transpose().mult(soln.F).mult(soln.C)
				.plus(soln.Ct.mult(Fderiv).mult(soln.C)).plus(soln.Ct.mult(soln.F).mult(Caderiv));

		return Eaderiv.diag().get(soln.getRm().nOccAlpha - 1);
	}

	public static double MNDOHomoDerivNew(SolutionU soln, SimpleMatrix xa, SimpleMatrix Faderiv) {

		SimpleMatrix Caderiv = soln.ca.transpose().mult(xa);
		SimpleMatrix Eaderiv = Caderiv.transpose().mult(soln.Fa).mult(soln.ca.transpose())
				.plus(soln.ca.mult(Faderiv).mult(soln.ca.transpose())).plus(soln.ca.mult(soln.Fa).mult(Caderiv));

		return Eaderiv.diag().get(soln.getRm().nOccAlpha - 1);
	}

	public static double MNDODipoleDeriv(Solution soln,
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
					D1deriv = atoms[i].D1pd(paramnum - 5);
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
			com[0] = com[0] + atom.getAtomProperties().getMass() * atom.getCoordinates()[0];
			com[1] = com[1] + atom.getAtomProperties().getMass() * atom.getCoordinates()[1];
			com[2] = com[2] + atom.getAtomProperties().getMass() * atom.getCoordinates()[2];
			mass += atom.getAtomProperties().getMass();
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
				hybridip[0] = hybridip[0] - 2.5416 * 2 * atoms[j].D1() *
						densityderiv.get(index[j][0], index[j][1]);
				hybridip[1] = hybridip[1] - 2.5416 * 2 * atoms[j].D1() *
						densityderiv.get(index[j][0], index[j][2]);
				hybridip[2] = hybridip[2] - 2.5416 * 2 * atoms[j].D1() *
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
}
