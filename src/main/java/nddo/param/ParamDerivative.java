package nddo.param;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;

import static nddo.State.nom;

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

		throw new IllegalArgumentException("Model is not MNDO!");
	}

	private static double alphaHfderiv(Solution soln, int Z) {
		double sum = 0;

		for (int i = 0; i < soln.atoms.length; i++) {
			for (int j = i + 1; j < soln.atoms.length; j++) {
				sum += soln.atoms[i].crfalphapd(soln.atoms[j],
						getNum(soln.atomicNumbers[i], soln.atomicNumbers[j], Z));
			}
		}

		return sum / Constants.HEATCONV;
	}

	private static double betaHfderiv(Solution soln, int Z, int type) {
		SimpleMatrix densitymatrix = soln.densityMatrix();

		NDDOOrbital[] orbitals = soln.orbitals;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomOfOrb[j] != atomOfOrb[k]) {
					double H = nom.Hbetapd(orbitals[j], orbitals[k],
							getNumBeta(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]],
									Z, orbitals[j].getL(), orbitals[k].getL(), type));
					e += 2 * densitymatrix.get(j, k) * H;
				}
			}
		}

		return e / Constants.HEATCONV;
	}

	private static double uxxHfderiv(Solution soln, int Z, int type) {
		SimpleMatrix densitymatrix = soln.densityMatrix();

		NDDOOrbital[] orbitals = soln.orbitals;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

		double e = 0;
		for (int j = 0; j < orbitals.length; j++) {
			if (atomicNumbers[atomOfOrb[j]] == Z && orbitals[j].getL() == type) {
				e += densitymatrix.get(j, j);
			}
		}

		return e / Constants.HEATCONV;
	}

	public static double zetaHfderiv(SolutionR soln, int Z, int type) {
		SimpleMatrix densitymatrix = soln.densityMatrix();
		NDDOOrbital[] orbitals = soln.orbitals;

		int[][] orbsOfAtom = soln.orbsOfAtom;
		int[][] missingIndex = soln.missingOfAtom;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

		SimpleMatrix H = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				if (atomOfOrb[j] == atomOfOrb[k]) {
					double Huv = 0;
					for (int an = 0; an < soln.atoms.length; an++)
						if (atomOfOrb[j] != an) Huv += soln.atoms[an].Vpd(orbitals[j], orbitals[k],
								getNum(atomicNumbers[an], atomicNumbers[atomOfOrb[j]], Z), type);
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { /* case 3*/
					double Huk = nom.Hzetapd(orbitals[j], orbitals[k],
							getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z), type);
					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}

		SimpleMatrix G = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double sum = 0;

				if (atomOfOrb[j] == atomOfOrb[k]) {
					for (int l : missingIndex[atomOfOrb[j]]) {
						for (int m : missingIndex[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								sum += soln.densityMatrix().get(l, m) *
										nom.Gpd(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
												getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[l]], Z),
												type);
							}
						}
					}
				}
				else for (int l : orbsOfAtom[atomOfOrb[j]]) {
					for (int m : orbsOfAtom[atomOfOrb[k]]) {
						sum += soln.densityMatrix().get(l, m) * (-0.5 *
								nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
										getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z), type));
					}
				}

				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		}

		return densitymatrix.elementMult(G.plusi(2, H)).elementSum() * 0.5 / Constants.HEATCONV;
	}

	private static double zetaHfderiv(SolutionU soln, int Z, int type) {
		NDDOOrbital[] orbitals = soln.orbitals;
		int[][] orbsOfAtom = soln.orbsOfAtom;
		int[][] missingOfAtom = soln.missingOfAtom;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

		SimpleMatrix H = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				if (atomOfOrb[j] == atomOfOrb[k]) {
					double Huv = 0;

					for (int an = 0; an < soln.atoms.length; an++) {
						if (atomOfOrb[j] != an) {
							Huv += soln.atoms[an].Vpd(orbitals[j], orbitals[k],
									getNum(atomicNumbers[an], atomicNumbers[atomOfOrb[j]], Z), type);
						}
					}
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { // case 3
					double Huk = nom.Hzetapd(orbitals[j], orbitals[k],
							getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z), type);
					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}

		SimpleMatrix Ga = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		SimpleMatrix Gb = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double suma = 0;
				double sumb = 0;

				if (atomOfOrb[j] == atomOfOrb[k]) {
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								suma += soln.densityMatrix().get(l, m) *
										nom.Gpd(orbitals[j], orbitals[k], orbitals[l],
												orbitals[m], getNum(atomicNumbers[atomOfOrb[j]],
														atomicNumbers[atomOfOrb[l]], Z), type);
								sumb += soln.densityMatrix().get(l, m) *
										nom.Gpd(orbitals[j], orbitals[k], orbitals[l],
												orbitals[m], getNum(atomicNumbers[atomOfOrb[j]],
														atomicNumbers[atomOfOrb[l]], Z), type);

							}
						}
					}
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							suma -= soln.alphaDensity().get(l, m) *
									nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
											getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z),
											type);

							sumb -= soln.betaDensity().get(l, m) *
									nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
											getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z),
											type);

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

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = 0; k < soln.nOrbitals; k++) {
				e += 0.5 * soln.alphaDensity().get(j, k) * (H.get(j, k) + Fa.get(j, k));
				e += 0.5 * soln.betaDensity().get(j, k) * (H.get(j, k) + Fb.get(j, k));

			}
		}

		double e2 = soln.alphaDensity().elementMult(Ga.plusi(H)).elementSum() +
				soln.betaDensity().elementMult(Gb.plusi(H)).elementSum() * 0.5 / Constants.HEATCONV;
		System.out.println("e2 = " + e2);

		return e / Constants.HEATCONV;

	}

	private static double eisolHfderiv(NDDOAtom[] atoms, int Z) {
		int counter = 0;

		for (NDDOAtom a : atoms) {
			if (a.getAtomProperties().getZ() == Z) {
				counter++;
			}
		}

		return -counter / Constants.HEATCONV;
	}

	public static SimpleMatrix[][] MNDOStaticMatrixDeriv(SolutionR soln, int Z, int firstParamIndex) {
		NDDOAtom[] atoms = soln.atoms;
		SimpleMatrix[] HDerivs = new SimpleMatrix[8];
		SimpleMatrix[] FDerivs = new SimpleMatrix[8];

		if (firstParamIndex <= 1) HDerivs[1] = betafockderivstatic(soln, Z, 0);
		if (firstParamIndex <= 1) FDerivs[1] = HDerivs[1].copy();
		if (firstParamIndex <= 3) HDerivs[3] = uxxfockderivstatic(soln, Z, 0);
		if (firstParamIndex <= 3) FDerivs[3] = HDerivs[3].copy();
		if (firstParamIndex <= 5) HDerivs[5] = zetaHderivstatic(atoms, soln, Z, 0);
		if (firstParamIndex <= 5) FDerivs[5] = HDerivs[5].copy().plus(zetaGderivstatic(atoms, soln, Z, 0));

		if (Z != 1) {
			if (firstParamIndex <= 2) HDerivs[2] = betafockderivstatic(soln, Z, 1);
			if (firstParamIndex <= 2) FDerivs[2] = HDerivs[2].copy();
			if (firstParamIndex <= 4) HDerivs[4] = uxxfockderivstatic(soln, Z, 1);
			if (firstParamIndex <= 4) FDerivs[4] = HDerivs[4].copy();
			if (firstParamIndex <= 6) HDerivs[6] = zetaHderivstatic(atoms, soln, Z, 1);
			if (firstParamIndex <= 6) FDerivs[6] = HDerivs[6].copy().plus(zetaGderivstatic(atoms, soln, Z, 1));
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

		return e / Constants.HEATCONV;
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
						sum -= (soln.Cta.get(i, u) * soln.Cta.get(j + NOccAlpha, v) +
								soln.Cta.get(j + NOccAlpha, u) * soln.Cta.get(i, v)) *
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
						sum -= (soln.Ctb.get(i, u) * soln.Ctb.get(j + NOccBeta, v) +
								soln.Ctb.get(j + NOccBeta, u) * soln.Ctb.get(i, v)) *
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

		SimpleMatrix Caderiv = soln.Cta.transpose().mult(xa);
		SimpleMatrix Eaderiv = Caderiv.transpose().mult(soln.Fa).mult(soln.Cta.transpose())
				.plus(soln.Cta.mult(Faderiv).mult(soln.Cta.transpose())).plus(soln.Cta.mult(soln.Fa).mult(Caderiv));

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
