package nddo.param;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import static nddo.State.nom;

public class ParamDerivative {
	public static double HfDeriv(Solution soln, int Z, int paramnum) {
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
			return eisolHfderiv(soln, Z);
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

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
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
		for (int j = 0; j < soln.nOrbitals; j++) {
			if (atomicNumbers[atomOfOrb[j]] == Z && orbitals[j].getL() == type) {
				e += densitymatrix.get(j, j);
			}
		}

		return e / Constants.HEATCONV;
	}

	private static double zetaHfderiv(SolutionR soln, int Z, int type) {
		SimpleMatrix densityMatrix = soln.densityMatrix();
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

					for (int an = 0; an < soln.atoms.length; an++) {
						if (atomOfOrb[j] != an) {
							Huv += soln.atoms[an].Vpd(orbitals[j], orbitals[k],
									getNum(atomicNumbers[an], atomicNumbers[atomOfOrb[j]], Z), type);
						}
					}

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

		return G.plusi(2, H).elementMulti(densityMatrix).elementSum() * 0.5 / Constants.HEATCONV;
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

		SimpleMatrix densityMatrix = soln.densityMatrix();
		SimpleMatrix alphaDensity = soln.alphaDensity();
		SimpleMatrix betaDensity = soln.betaDensity();

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double suma = 0;
				double sumb = 0;

				if (atomOfOrb[j] == atomOfOrb[k]) {
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								double v = densityMatrix.get(l, m) *
										nom.Gpd(orbitals[j], orbitals[k], orbitals[l],
												orbitals[m], getNum(atomicNumbers[atomOfOrb[j]],
														atomicNumbers[atomOfOrb[l]], Z), type);

								suma += v;
								sumb += v;
							}
						}
					}
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							double gpd = nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
									getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z),
									type);

							suma -= alphaDensity.get(l, m) * gpd;
							sumb -= betaDensity.get(l, m) * gpd;
						}
					}
				}

				Ga.set(j, k, suma);
				Ga.set(k, j, suma);

				Gb.set(j, k, sumb);
				Gb.set(k, j, sumb);
			}
		}

		H.scalei(2);

		return (Ga.plusi(H).elementMulti(alphaDensity).elementSum() +
				Gb.plusi(H).elementMulti(betaDensity).elementSum()) * 0.5 / Constants.HEATCONV;
	}

	private static double eisolHfderiv(Solution soln, int Z) {
		int counter = 0;

		for (NDDOAtom a : soln.atoms) {
			if (a.getAtomProperties().getZ() == Z) {
				counter++;
			}
		}

		return -counter / Constants.HEATCONV;
	}


	public static SimpleMatrix[][] staticDeriv(SolutionR soln, int Z, int firstParamIndex) {
		SimpleMatrix[] HDerivs = new SimpleMatrix[8];
		SimpleMatrix[] FDerivs = new SimpleMatrix[8];

		if (firstParamIndex <= 1) {
			FDerivs[1] = HDerivs[1] = betafockderivstatic(soln, Z, 0);
		}
		if (firstParamIndex <= 3) {
			FDerivs[3] = HDerivs[3] = uxxfockderivstatic(soln, Z, 0);
		}
		if (firstParamIndex <= 5) {
			HDerivs[5] = zetaHderivstatic(soln, Z, 0);
			FDerivs[5] = HDerivs[5].plus(zetaGderivstatic(soln, Z, 0));
		}

		if (Z != 1) {
			if (firstParamIndex <= 2) {
				FDerivs[2] = HDerivs[2] = betafockderivstatic(soln, Z, 1);
			}
			if (firstParamIndex <= 4) {
				FDerivs[4] = HDerivs[4] = uxxfockderivstatic(soln, Z, 1);
			}
			if (firstParamIndex <= 6) {
				HDerivs[6] = zetaHderivstatic(soln, Z, 1);
				FDerivs[6] = HDerivs[6].plus(zetaGderivstatic(soln, Z, 1));
			}
		}

		return new SimpleMatrix[][]{HDerivs, FDerivs};
	}

	public static SimpleMatrix[][] staticDeriv(SolutionU soln, int Z, int firstParamIndex) {
		SimpleMatrix[] HDerivs = new SimpleMatrix[8];
		SimpleMatrix[] FaDerivs = new SimpleMatrix[8];
		SimpleMatrix[] FbDerivs = new SimpleMatrix[8];

		if (firstParamIndex <= 1) {
			FaDerivs[1] = FbDerivs[1] = HDerivs[1] = betafockderivstatic(soln, Z, 0);
		}
		if (firstParamIndex <= 3) {
			FbDerivs[3] = FaDerivs[3] = HDerivs[3] = uxxfockderivstatic(soln, Z, 0);
		}
		if (firstParamIndex <= 5) {
			HDerivs[5] = zetaHderivstatic(soln, Z, 0);

			SimpleMatrix[] matrices = zetaGderivstatic(soln, Z, 0);
			FaDerivs[5] = HDerivs[5].plus(matrices[0]);
			FbDerivs[5] = HDerivs[5].plus(matrices[1]);
		}

		if (Z != 1) {
			if (firstParamIndex <= 2) {
				FbDerivs[2] = FaDerivs[2] = HDerivs[2] = betafockderivstatic(soln, Z, 1);
			}
			if (firstParamIndex <= 4) {
				FbDerivs[4] = FaDerivs[4] = HDerivs[4] = uxxfockderivstatic(soln, Z, 1);
			}
			if (firstParamIndex <= 6) {
				HDerivs[6] = zetaHderivstatic(soln, Z, 1);

				SimpleMatrix[] matrices = zetaGderivstatic(soln, Z, 1);
				FaDerivs[6] = HDerivs[6].plus(matrices[0]);
				FbDerivs[6] = HDerivs[6].plus(matrices[1]);
			}
		}

		return new SimpleMatrix[][]{HDerivs, FaDerivs, FbDerivs};
	}

	private static SimpleMatrix betafockderivstatic(Solution soln, int Z, int type) {
		SimpleMatrix F = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		NDDOOrbital[] orbitals = soln.orbitals;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				if (atomOfOrb[j] != atomOfOrb[k]) {
					double H = nom.Hbetapd(orbitals[j], orbitals[k],
							getNumBeta(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]],
									Z, orbitals[j].getL(), orbitals[k].getL(), type));
					F.set(j, k, H);
					F.set(k, j, H);
				}
			}
		}

		return F;
	}

	private static SimpleMatrix uxxfockderivstatic(Solution soln, int Z, int type) {
		SimpleMatrix F = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int j = 0; j < soln.nOrbitals; j++) {
			if (soln.atomicNumbers[soln.atomOfOrb[j]] == Z && soln.orbitals[j].getL() == type) {
				F.set(j, j, 1);
			}
		}

		return F;
	}

	private static SimpleMatrix zetaHderivstatic(Solution soln, int Z, int type) {
		NDDOOrbital[] orbitals = soln.orbitals;
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

		return H;
	}

	private static SimpleMatrix zetaGderivstatic(SolutionR soln, int Z, int type) { // todo unify first and second
		// order
		NDDOOrbital[] orbitals = soln.orbitals;
		int[][] orbsOfAtom = soln.orbsOfAtom;
		int[][] missingOfAtom = soln.missingOfAtom;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

		SimpleMatrix G = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double sum = 0;

				if (atomOfOrb[j] == atomOfOrb[k]) {
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								sum += soln.densityMatrix().get(l, m) *
										nom.Gpd(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
												getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[l]], Z),
												type);

							}
						}
					}
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							sum += soln.densityMatrix().get(l, m) * (-0.5 *
									nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
											getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z),
											type));
						}
					}
				}

				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		}


		return G;
	}

	private static SimpleMatrix[] zetaGderivstatic(SolutionU soln, int Z, int type) {
		NDDOOrbital[] orbitals = soln.orbitals;
		int[][] orbsOfAtom = soln.orbsOfAtom;
		int[][] missingOfAtom = soln.missingOfAtom;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

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
										State.nom.Gpd(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
												getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[l]], Z),
												type);

								sumb += soln.densityMatrix().get(l, m) *
										State.nom.Gpd(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
												getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[l]], Z),
												type);
							}
						}
					}
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							suma -= soln.alphaDensity().get(l, m) *
									State.nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
											getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z), type);

							sumb -= soln.betaDensity().get(l, m) *
									State.nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
											getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z), type);
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


	public static SimpleMatrix xMatrix(SolutionR soln, SimpleMatrix F) {
		int numRows = F.numRows();
		SimpleMatrix x = new SimpleMatrix(numRows, numRows);

		for (int i = 0; i < numRows; i++) {
			for (int j = i + 1; j < numRows; j++) {
				double ej = soln.E.get(j);
				double ei = soln.E.get(i);
				
				double v = ej == ei ? 0 : F.get(i, j) / (ej - ei);

				x.set(i, j, v);
				x.set(j, i, -v);
			}
		}

		return x;
	}

	public static SimpleMatrix[] xMatrix(SolutionU soln, SimpleMatrix Fa, SimpleMatrix Fb) {
		int numRows = Fa.numRows();
		SimpleMatrix xa = new SimpleMatrix(numRows, numRows);
		SimpleMatrix xb = new SimpleMatrix(numRows, numRows);

		for (int i = 0; i < numRows; i++) {
			for (int j = i + 1; j < numRows; j++) {
				double eaj = soln.Ea.get(j);
				double eai = soln.Ea.get(i);
				double ebj = soln.Eb.get(j);
				double ebi = soln.Eb.get(i);

				double va = eaj == eai ? 0 : Fa.get(i, j) / (eaj - eai);
				double vb = ebj == ebi ? 0 : Fb.get(i, j) / (ebj - ebi);

				xa.set(i, j, va);
				xa.set(j, i, -va);

				xb.set(i, j, vb);
				xb.set(j, i, -vb);
			}
		}

		return new SimpleMatrix[]{xa, xb};
	}


	/**
	 * Computes Hf faster if Hderiv and Fderiv are already available.
	 */
	public static double HfDeriv(SolutionR soln, SimpleMatrix Hderiv, SimpleMatrix Fderiv) {
		double e = 0;

		SimpleMatrix densitymatrix = soln.densityMatrix();

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = 0; k < soln.nOrbitals; k++) {
				// cannot be converted to SimpleMatrix operations whilst preserving H and Fderiv
				e += densitymatrix.get(j, k) * (Hderiv.get(j, k) + Fderiv.get(j, k));
			}
		}

		return e * 0.5 / Constants.HEATCONV;
	}

	public static double HfDeriv(SolutionU soln, SimpleMatrix Hderiv, SimpleMatrix Faderiv, SimpleMatrix Fbderiv) {
		double e = 0;

		SimpleMatrix alphaDensity = soln.alphaDensity();
		SimpleMatrix betaDensity = soln.betaDensity();

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = 0; k < soln.nOrbitals; k++) {
				double v = Hderiv.get(j, k);
				e += alphaDensity.get(j, k) * (v + Faderiv.get(j, k)) +
						betaDensity.get(j, k) * (v + Fbderiv.get(j, k));
			}
		}

		return e * 0.5 / Constants.HEATCONV;
	}

	public static double nddoDipoleDeriv(Solution soln, SimpleMatrix densityDeriv, int Z, int paramnum) {
		if (soln.dipole <= 1E-6) {
			return 0;
		}

		NDDOAtom[] atoms = soln.atoms;
		int[][] orbsOfAtom = soln.orbsOfAtom;
		SimpleMatrix densityMatrix = soln.densityMatrix();
		double D1deriv = 0;

		if (paramnum == 5 || paramnum == 6) {
			for (NDDOAtom atom : atoms) {
				if (atom.getAtomProperties().getZ() == Z) {
					D1deriv = atom.D1pd(paramnum - 5);
					break;
				}
			}
		}

		double[] populationsderiv = new double[atoms.length];
		for (int j = 0; j < atoms.length; j++) {
			double sum = 0;

			for (int k : orbsOfAtom[j]) {
				sum += densityDeriv.get(k, k);
			}

			populationsderiv[j] = -sum;
		}


		double[] com = new double[]{0, 0, 0};
		double mass = 0;
		for (NDDOAtom atom : atoms) {
			com[0] += atom.getAtomProperties().getMass() * atom.getCoordinates()[0];
			com[1] += atom.getAtomProperties().getMass() * atom.getCoordinates()[1];
			com[2] += atom.getAtomProperties().getMass() * atom.getCoordinates()[2];

			mass += atom.getAtomProperties().getMass();
		}
		com[0] /= mass;
		com[1] /= mass;
		com[2] /= mass;

		double[] chargedip = new double[]{0, 0, 0};
		for (int j = 0; j < atoms.length; j++) {
			chargedip[0] += Constants.DIPOLECONV * populationsderiv[j] * (atoms[j].getCoordinates()[0] - com[0]);
			chargedip[1] += Constants.DIPOLECONV * populationsderiv[j] * (atoms[j].getCoordinates()[1] - com[1]);
			chargedip[2] += Constants.DIPOLECONV * populationsderiv[j] * (atoms[j].getCoordinates()[2] - com[2]);
		}

		double[] hybridip = new double[]{0, 0, 0};
		for (int j = 0; j < atoms.length; j++) {
			if (orbsOfAtom[j].length > 1) { // exclude hydrogen
				hybridip[0] -=
						Constants.DIPOLECONV * 2 * atoms[j].D1() * densityDeriv.get(orbsOfAtom[j][0],
								orbsOfAtom[j][1]);
				hybridip[1] -=
						Constants.DIPOLECONV * 2 * atoms[j].D1() * densityDeriv.get(orbsOfAtom[j][0],
								orbsOfAtom[j][2]);
				hybridip[2] -=
						Constants.DIPOLECONV * 2 * atoms[j].D1() * densityDeriv.get(orbsOfAtom[j][0],
								orbsOfAtom[j][3]);

				if (atoms[j].getAtomProperties().getZ() == Z) {
					hybridip[0] -=
							Constants.DIPOLECONV * 2 * D1deriv * densityMatrix.get(orbsOfAtom[j][0], orbsOfAtom[j][1]);
					hybridip[1] -=
							Constants.DIPOLECONV * 2 * D1deriv * densityMatrix.get(orbsOfAtom[j][0], orbsOfAtom[j][2]);
					hybridip[2] -=
							Constants.DIPOLECONV * 2 * D1deriv * densityMatrix.get(orbsOfAtom[j][0], orbsOfAtom[j][3]);
				}
			}
		}

		double[] dipoletot =
				new double[]{chargedip[0] + hybridip[0], chargedip[1] + hybridip[1], chargedip[2] + hybridip[2]};

		return (dipoletot[0] * soln.dipoletot[0] +
				dipoletot[1] * soln.dipoletot[1] +
				dipoletot[2] * soln.dipoletot[2]) / soln.dipole;
	}

	public static double homoDeriv(SolutionR soln, SimpleMatrix x, SimpleMatrix Fderivmo) {
		SimpleMatrix Cderiv = soln.C.mult(x); // todo are you sure this is faster?
		SimpleMatrix Ederiv = Utils.plusTrans(soln.Ct.mult(soln.F).mult(Cderiv)).plusi(Fderivmo);

		return Ederiv.get(soln.rm.nOccAlpha - 1, soln.rm.nOccAlpha - 1);
	}

	public static double homoDeriv(SolutionU soln, SimpleMatrix xa, SimpleMatrix Faderivmo) {
		SimpleMatrix Cderiv = soln.Ca.mult(xa);
		SimpleMatrix Ederiv = Utils.plusTrans(soln.Cta.mult(soln.Fa).mult(Cderiv)).plusi(Faderivmo);

		return Ederiv.get(soln.rm.nOccAlpha - 1, soln.rm.nOccAlpha - 1);
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


	@Deprecated
	public static double MNDOHomoDerivtemp(SolutionR soln, SimpleMatrix coeffDeriv, SimpleMatrix Fderiv) {

		int index = (int) (soln.nElectrons / 2.0) - 1;

		if (index < 0) {
			return 0;
		}

		SimpleMatrix coeff = soln.Ct.extractVector(true, index);

		double sum = 0;

		for (int i = 0; i < soln.nOrbitals; i++) {
			for (int j = 0; j < soln.nOrbitals; j++) {
				sum += soln.F.get(i, j) * (coeff.get(i) * coeffDeriv.get(j) +
						coeff.get(j) * coeffDeriv.get(i));
				sum += Fderiv.get(i, j) * coeff.get(i) * coeff.get(j);
			}
		}

		return sum;


	}

	@Deprecated
	public static SimpleMatrix xArrayComplementary(SolutionR soln, SimpleMatrix fockderiv) {
		int NOcc = (int) (soln.nElectrons / 2.0);
		if (NOcc == 0) {
			return new SimpleMatrix(0, 0);
		}
		SimpleMatrix x = new SimpleMatrix(NOcc - 1, 1);
		int count1 = 0;

		for (int j = 0; j < NOcc - 1; j++) {
			double element = 0;

			for (int u = 0; u < soln.nOrbitals; u++) {
				for (int v = 0; v < soln.nOrbitals; v++) {
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

	@Deprecated
	public static SimpleMatrix xarrayForIE(SolutionR soln, SimpleMatrix xlimited, SimpleMatrix xcomplementary) {

		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.nOrbitals - NOcc;

		SimpleMatrix x = new SimpleMatrix(soln.nOrbitals - 1, 1);

		int count = 0;

		for (int i = xlimited.getNumElements() - NVirt; i < xlimited.getNumElements(); i++) {
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

	@Deprecated
	public static SimpleMatrix HOMOCoefficientDerivativeComplementary(SimpleMatrix x, SolutionR soln) {


		SimpleMatrix CDeriv = new SimpleMatrix(1, soln.nOrbitals);

		int NOcc = (int) (soln.nElectrons / 2.0);

		for (int u = 0; u < soln.nOrbitals; u++) {
			double sum = 0;

			for (int k = 0; k < soln.nOrbitals; k++) {

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
}
