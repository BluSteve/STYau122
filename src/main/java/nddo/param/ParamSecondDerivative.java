package nddo.param;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.geometry.GeometryDerivative;
import nddo.math.PopleThiel;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.apache.commons.lang3.time.StopWatch;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import tools.Pow;
import tools.Utils;

import java.util.ArrayList;

public class ParamSecondDerivative {
	public static double alphaHfderiv2(Solution soln, int Z1) {
		double sum = 0;

		for (int i = 0; i < soln.atoms.length; i++) {
			for (int j = i + 1; j < soln.atoms.length; j++) {
				sum += soln.atoms[i].crfalphap2d(soln.atoms[j],
						ParamDerivative.getNum(soln.atomicNumbers[i], soln.atomicNumbers[j], Z1));
			}
		}

		return sum / Constants.HEATCONV;
	}

	public static SimpleMatrix gammaMatrix(SolutionR soln, SimpleMatrix totalderiv) {
		int numRows = totalderiv.numRows();
		SimpleMatrix gamma = new SimpleMatrix(numRows, numRows);

		for (int i = 0; i < numRows; i++) {
			for (int j = i + 1; j < numRows; j++) {
				double ej = soln.E.get(j);
				double ei = soln.E.get(i);

				double va = ej == ei ? 0 : totalderiv.get(i, j) / (ej - ei);
				double vb = ej == ei ? 0 : totalderiv.get(j, i) / (ei - ej);

				gamma.set(i, j, va);
				gamma.set(j, i, vb);
			}
		}

		return gamma;

	}

	public static SimpleMatrix gammaMatrix(SolutionU soln, SimpleMatrix totalderivalpha) {
		int numRows = totalderivalpha.numRows();
		SimpleMatrix gamma = new SimpleMatrix(numRows, numRows);

		for (int i = 0; i < numRows; i++) {
			for (int j = i + 1; j < numRows; j++) {
				double ej = soln.Ea.get(j);
				double ei = soln.Ea.get(i);

				double va = ej == ei ? 0 : totalderivalpha.get(i, j) / (ej - ei);
				double vb = ej == ei ? 0 : totalderivalpha.get(j, i) / (ei - ej);

				gamma.set(i, j, va);
				gamma.set(j, i, vb);
			}
		}

		return gamma;
	}


	private static SimpleMatrix densityDeriv2static(int nOrbitals, int nOcc, SimpleMatrix C, SimpleMatrix Ct,
													SimpleMatrix COcc, SimpleMatrix CtOcc, SimpleMatrix x1,
													SimpleMatrix x2) {
		SimpleMatrix x1occ = x1.extractMatrix(0, nOrbitals, 0, nOcc);
		SimpleMatrix x2occ = x2.extractMatrix(0, nOrbitals, 0, nOcc);

		SimpleMatrix mat = Utils.plusTrans(C.mult(x2.mult(x1occ).plusi(x1.mult(x2occ))).mult(CtOcc));
		SimpleMatrix mat2 = Utils.plusTrans(C.mult(x1occ).mult(x2occ.transpose()).mult(Ct));
		SimpleMatrix mat3 = Utils.plusTrans(COcc.mult(x1occ.transpose()).mult(x2occ).mult(CtOcc));

		return mat.plusi(mat2).plusi(mat3);
	}

	public static SimpleMatrix densityDeriv2static(SolutionR s, SimpleMatrix x1, SimpleMatrix x2) {
		return densityDeriv2static(s.nOrbitals, s.rm.nOccAlpha, s.C, s.Ct, s.COcc, s.CtOcc, x1, x2).scalei(2);
	}

	public static SimpleMatrix[] densityDeriv2static(SolutionU s, SimpleMatrix[] x1, SimpleMatrix[] x2) {
		SimpleMatrix sm =
				densityDeriv2static(s.nOrbitals, s.rm.nOccAlpha, s.Ca, s.Cta, s.CaOcc, s.CtaOcc, x1[0], x2[0]);
		SimpleMatrix sm2 =
				densityDeriv2static(s.nOrbitals, s.rm.nOccBeta, s.Cb, s.Ctb, s.CbOcc, s.CtbOcc, x1[1], x2[1]);

		return new SimpleMatrix[]{sm, sm2};
	}


	public static SimpleMatrix Gderivstatic(SolutionR soln, SimpleMatrix mat, int Z, int paramNum) {
		if (paramNum < 5 || paramNum > 6) {
			return new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		}
		else {
			return zetaGderivstatic(soln, mat, Z, paramNum - 5);
		}
	}

	public static SimpleMatrix[] Gderivstatic(SolutionU soln, SimpleMatrix alphamat, SimpleMatrix betamat, int Z,
											  int paramNum) {
		if (paramNum < 5 || paramNum > 6) {
			return new SimpleMatrix[]{new SimpleMatrix(soln.nOrbitals, soln.nOrbitals),
					new SimpleMatrix(soln.nOrbitals, soln.nOrbitals)};
		}
		else {
			return zetaGderivstatic(soln, alphamat, betamat, Z, paramNum - 5);
		}
	}

	private static SimpleMatrix zetaGderivstatic(SolutionR soln, SimpleMatrix mat, int Z, int type) {
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
								sum += mat.get(l, m) * State.nom.Gpd(orbitals[j], orbitals[k], orbitals[l],
										orbitals[m],
										ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]],
												atomicNumbers[atomOfOrb[l]],
												Z), type);
							}
						}
					}
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							sum += mat.get(l, m) * (-0.5 *
									State.nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
											ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]],
													atomicNumbers[atomOfOrb[k]], Z), type));
						}
					}
				}

				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		}

		return G;
	}

	private static SimpleMatrix[] zetaGderivstatic(SolutionU soln, SimpleMatrix alphamat, SimpleMatrix betamat, int Z,
												   int type) {
		NDDOOrbital[] orbitals = soln.orbitals;
		int[][] orbsOfAtom = soln.orbsOfAtom;
		int[][] missingOfAtom = soln.missingOfAtom;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

		SimpleMatrix Ga = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		SimpleMatrix Gb = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		SimpleMatrix mat = alphamat.plus(betamat);

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double suma = 0;
				double sumb = 0;

				if (atomOfOrb[j] == atomOfOrb[k]) {
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								double v = mat.get(l, m) *
										State.nom.Gpd(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
												ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]],
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
							double v = State.nom.Gpd(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
									ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]],
											Z),
									type);

							suma -= alphamat.get(l, m) * v;
							sumb -= betamat.get(l, m) * v;
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


	public static SimpleMatrix Hderiv2(Solution soln, int Z1, int param1, int Z2, int param2) {
		if (param1 >= 1 && param1 <= 2) {//beta
			if (param2 >= 5 && param2 <= 6) {//zeta
				return ParamSecondDerivative.betazetaHderiv2(soln, Z1, param1 - 1, Z2, param2 - 5);
			}
		}

		else if (param1 >= 5 && param1 <= 6) {//zeta
			if (param2 >= 1 && param2 <= 2) {//beta
				return ParamSecondDerivative.betazetaHderiv2(soln, Z2, param2 - 1, Z1, param1 - 5);
			}
			else if (param2 >= 5 && param2 <= 6) {
				return ParamSecondDerivative.zetazetaHderiv2(soln, Z1, param1 - 5, Z2, param2 - 5);
			}
		}

		return new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
	}

	protected static SimpleMatrix betazetaHderiv2(Solution soln, int Z1, int type1, int Z2, int type2) {
		NDDOOrbital[] orbitals = soln.orbitals;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

		SimpleMatrix H = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				if (atomOfOrb[j] != atomOfOrb[k]) {
					double Huk = State.nom.Hbetazetap2d(orbitals[j], orbitals[k],
							ParamDerivative.getNumBeta(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z1,
									orbitals[j].getL(), orbitals[k].getL(), type1),
							ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z2),
							type2);

					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}

		return H;
	}

	protected static SimpleMatrix zetazetaHderiv2(Solution soln, int Z1, int type1, int Z2, int type2) {
		NDDOAtom[] atoms = soln.atoms;
		NDDOOrbital[] orbitals = soln.orbitals;
		int[] atomOfOrb = soln.atomOfOrb;
		int[] atomicNumbers = soln.atomicNumbers;

		SimpleMatrix H = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				if (atomOfOrb[j] == atomOfOrb[k]) {
					double Huv = 0;

					for (int an = 0; an < atoms.length; an++) {
						if (atomOfOrb[j] != an) {
							Huv += atoms[an].Vp2d(orbitals[j], orbitals[k],
									ParamDerivative.getNum(atomicNumbers[an], atomicNumbers[atomOfOrb[j]], Z1), type1,
									ParamDerivative.getNum(atomicNumbers[an], atomicNumbers[atomOfOrb[j]], Z2), type2);
						}
					}

					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else {
					double Huk = State.nom.Hzetazetap2d(orbitals[j], orbitals[k],
							ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z1),
							type1,
							ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]], Z2),
							type2);

					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}

		return H;
	}


	public static SimpleMatrix Gderiv2static(SolutionR soln, int Z1, int param1, int Z2, int param2) {
		if (param1 >= 5 && param1 <= 6) { // zeta
			if (param2 >= 5 && param2 <= 6) {
				return ParamSecondDerivative.zetazetaGderiv2static(soln, Z1, param1 - 5, Z2, param2 - 5);
			}
		}
		return new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
	}

	public static SimpleMatrix[] Gderiv2static(SolutionU soln, int Z1, int param1, int Z2, int param2) {
		if (param1 >= 5 && param1 <= 6) { // zeta
			if (param2 >= 5 && param2 <= 6) {
				return ParamSecondDerivative.zetazetaGderiv2static(soln, Z1, param1 - 5, Z2, param2 - 5);
			}
		}
		return new SimpleMatrix[]{new SimpleMatrix(soln.nOrbitals, soln.nOrbitals),
				new SimpleMatrix(soln.nOrbitals, soln.nOrbitals)};
	}

	private static SimpleMatrix zetazetaGderiv2static(SolutionR soln, int Z1, int type1, int Z2, int type2) {
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
										State.nom.Gp2d(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
												ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]],
														atomicNumbers[atomOfOrb[l]], Z1), type1,
												ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]],
														atomicNumbers[atomOfOrb[l]], Z2), type2);

							}
						}
					}
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							sum += soln.densityMatrix().get(l, m) * (-0.5 *
									State.nom.Gp2d(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
											ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]],
													atomicNumbers[atomOfOrb[k]], Z1), type1,
											ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]],
													atomicNumbers[atomOfOrb[k]], Z2), type2));
						}
					}
				}

				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		}

		return G;
	}

	private static SimpleMatrix[] zetazetaGderiv2static(SolutionU soln, int Z1, int type1, int Z2, int type2) {
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
								double v = soln.densityMatrix().get(l, m) *
										State.nom.Gp2d(orbitals[j], orbitals[k], orbitals[l], orbitals[m],
												ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]],
														atomicNumbers[atomOfOrb[l]], Z1), type1,
												ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]],
														atomicNumbers[atomOfOrb[l]], Z2), type2);

								suma += v;
								sumb += v;

							}
						}
					}
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							double v = -State.nom.Gp2d(orbitals[j], orbitals[l], orbitals[k], orbitals[m],
									ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]],
											Z1), type1,
									ParamDerivative.getNum(atomicNumbers[atomOfOrb[j]], atomicNumbers[atomOfOrb[k]],
											Z2), type2);

							suma += soln.alphaDensity().get(l, m) * v;
							sumb += soln.betaDensity().get(l, m) * v;

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


	public static SimpleMatrix staticFockDeriv(SolutionR soln, SimpleMatrix Fstatictotal, SimpleMatrix densityA,
											   SimpleMatrix densityB, SimpleMatrix densityderiv2static, int Z1,
											   int param1, int Z2, int param2) {
		SimpleMatrix GaB = ParamSecondDerivative.Gderivstatic(soln, densityB, Z1, param1);
		SimpleMatrix GbA = ParamSecondDerivative.Gderivstatic(soln, densityA, Z2, param2);

		SimpleMatrix omega = PopleThiel.responseMatrix(soln, densityderiv2static);

		return omega.plusi(GaB).plusi(GbA).plusi(Fstatictotal);
	}

	public static SimpleMatrix[] staticFockDeriv(SolutionU soln, SimpleMatrix[] Fstatictotal,
												 SimpleMatrix[] densitiesA,
												 SimpleMatrix[] densitiesB, SimpleMatrix[] densityderivs2static,
												 int Z1,
												 int param1, int Z2, int param2) {
		SimpleMatrix[] GaB = ParamSecondDerivative.Gderivstatic(soln, densitiesB[0], densitiesB[1], Z1, param1);
		SimpleMatrix[] GbA = ParamSecondDerivative.Gderivstatic(soln, densitiesA[0], densitiesA[1], Z2, param2);
		SimpleMatrix[] omega = PopleThiel.responseMatrix(soln, densityderivs2static);

		SimpleMatrix Phialpha = omega[0].plusi(GaB[0]).plusi(GbA[0]).plusi(Fstatictotal[0]);
		SimpleMatrix Phibeta = omega[1].plusi(GaB[1]).plusi(GbA[1]).plusi(Fstatictotal[1]);

		return new SimpleMatrix[]{Phialpha, Phibeta};
	}

	private static SimpleMatrix staticMatrix(SimpleMatrix Ct, SimpleMatrix C, SimpleMatrix E, SimpleMatrix Phi,
											 SimpleMatrix FockA, SimpleMatrix FockB, SimpleMatrix xA,
											 SimpleMatrix xB) {
		SimpleMatrix FA = Ct.mult(FockA).mult(C);
		double[] arrfa = FA.diag().getDDRM().data;
		SimpleMatrix FB = Ct.mult(FockB).mult(C);
		double[] arrfb = FB.diag().getDDRM().data;
		SimpleMatrix matrix = Ct.mult(Phi).mult(C);

		SimpleMatrix xbFA = xB.copy();
		CommonOps_DDRM.multRows(arrfa, xbFA.getDDRM());
		Utils.plusTrans(xbFA);
		SimpleMatrix xaFB = xA.copy();
		CommonOps_DDRM.multRows(arrfb, xaFB.getDDRM());
		Utils.plusTrans(xaFB);
		matrix.plusi(xbFA).plusi(xaFB);

		double[] Earr = E.getDDRM().data;

		final SimpleMatrix xaE = xA.copy();
		CommonOps_DDRM.multRows(Earr, xaE.getDDRM());

		final SimpleMatrix xbE = xB.copy();
		CommonOps_DDRM.multRows(Earr, xbE.getDDRM());

		SimpleMatrix xatxb = xaE.transpose().mult(xB);
		Utils.plusTrans(xatxb);
		SimpleMatrix xaxb = xaE.mult(xB);
		SimpleMatrix xbxa = xbE.mult(xA);

		matrix.minusi(xatxb);
		matrix.minusi(xaxb).minusi(xbxa);

		return matrix;
	}

	public static SimpleMatrix staticMatrix(SolutionR soln, SimpleMatrix Phi, SimpleMatrix FockA, SimpleMatrix FockB,
											SimpleMatrix xA, SimpleMatrix xB) {
		return staticMatrix(soln.Ct, soln.C, soln.E, Phi, FockA, FockB, xA, xB);
	}

	public static SimpleMatrix[] staticMatrix(SolutionU soln, SimpleMatrix[] Phi, SimpleMatrix[] FockA,
											  SimpleMatrix[] FockB, SimpleMatrix[] xA, SimpleMatrix[] xB) {
		return new SimpleMatrix[]{staticMatrix(soln.Cta, soln.Ca, soln.Ea, Phi[0], FockA[0], FockB[0], xA[0], xB[0]),
				staticMatrix(soln.Ctb, soln.Cb, soln.Eb, Phi[1], FockA[1], FockB[1], xA[1], xB[1])};
	}


	public static double HfDeriv2(SolutionR soln, int Z1, int param1, int Z2, int param2, SimpleMatrix Hfirst,
								  SimpleMatrix Ffirst, SimpleMatrix densityderiv, int densityderivparamtype) {

		SimpleMatrix F = null;

		SimpleMatrix H = Hderiv2(soln, Z1, param1, Z2, param2);


		switch (densityderivparamtype) {
			case 0:
				F = H.plus(Gderiv2static(soln, Z1, param1, Z2, param2))
						.plus(Gderivstatic(soln, densityderiv, Z1, param1));
				break;
			case 1:
				F = H.plus(Gderiv2static(soln, Z1, param1, Z2, param2))
						.plus(Gderivstatic(soln, densityderiv, Z2, param2));
		}

		double e = 0;

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = 0; k < soln.nOrbitals; k++) {
				e += 0.5 * soln.densityMatrix().get(j, k) * (H.get(j, k) + F.get(j, k));
				e += 0.5 * densityderiv.get(j, k) * (Hfirst.get(j, k) + Ffirst.get(j, k));

			}
		}

		return e / Constants.HEATCONV;

	}

	public static double HfDeriv2(SolutionU soln, int Z1, int param1, int Z2, int param2, SimpleMatrix Hfirst,
								  SimpleMatrix Fafirst, SimpleMatrix Fbfirst, SimpleMatrix densityderivalpha,
								  SimpleMatrix densityderivbeta, int densityderivparamtype) {

		SimpleMatrix Fa = null;

		SimpleMatrix Fb = null;

		SimpleMatrix H = Hderiv2(soln, Z1, param1, Z2, param2);


		switch (densityderivparamtype) {
			case 0:

				SimpleMatrix[] matrices1 = Gderiv2static(soln, Z1, param1, Z2, param2);

				SimpleMatrix[] matrices2 = Gderivstatic(soln, densityderivalpha, densityderivbeta, Z1, param1);

				Fa = H.plus(matrices1[0]).plus(matrices2[0]);
				Fb = H.plus(matrices1[1]).plus(matrices2[1]);
				break;
			case 1:
				matrices1 = Gderiv2static(soln, Z1, param1, Z2, param2);

				matrices2 = Gderivstatic(soln, densityderivalpha, densityderivbeta, Z2, param2);
				Fa = H.plus(matrices1[0]).plus(matrices2[0]);
				Fb = H.plus(matrices1[1]).plus(matrices2[1]);
		}

		double e = 0;

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = 0; k < soln.nOrbitals; k++) {
				e += 0.5 * soln.alphaDensity().get(j, k) * (H.get(j, k) + Fa.get(j, k));
				e += 0.5 * densityderivalpha.get(j, k) * (Hfirst.get(j, k) + Fafirst.get(j, k));

				e += 0.5 * soln.betaDensity().get(j, k) * (H.get(j, k) + Fb.get(j, k));
				e += 0.5 * densityderivbeta.get(j, k) * (Hfirst.get(j, k) + Fbfirst.get(j, k));

			}
		}

		return e / Constants.HEATCONV;

	}

	public static double dipoleDeriv2(Solution soln, SimpleMatrix densityderiva, SimpleMatrix densityderivb,
									  SimpleMatrix densityderiv2, int Z1, int paramnum1, int Z2, int paramnum2) {

		if (soln.dipole <= 1E-6) {
			return 0;
		}

		NDDOAtom[] atoms = soln.atoms;

		double D1deriva = 0;

		double D1derivb = 0;

		double D1deriv2 = 0;

		if (paramnum1 == 5 || paramnum1 == 6) {
			for (int i = 0; i < atoms.length; i++) {
				if (atoms[i].getAtomProperties().getZ() == Z1) {
					D1deriva = atoms[i].D1pd(paramnum1 - 5);
					break;
				}
			}
		}

		if (paramnum2 == 5 || paramnum2 == 6) {
			for (int i = 0; i < atoms.length; i++) {
				if (atoms[i].getAtomProperties().getZ() == Z2) {
					D1derivb = atoms[i].D1pd(paramnum2 - 5);
					break;
				}
			}
		}

		if (Z1 == Z2 && (paramnum1 == 5 || paramnum1 == 6) && (paramnum2 == 5 || paramnum2 == 6)) {
			for (int i = 0; i < atoms.length; i++) {
				if (atoms[i].getAtomProperties().getZ() == Z1) {
					D1deriv2 = atoms[i].D1p2d(paramnum2 + paramnum1 - 10);
					break;
				}
			}
		}

		int[][] orbsOfAtom = soln.orbsOfAtom;

		SimpleMatrix densityMatrix = soln.densityMatrix();

		double[] populationsderiva = new double[atoms.length];
		double[] populationsderivb = new double[atoms.length];
		double[] populationsderiv2 = new double[atoms.length];


		for (int j = 0; j < atoms.length; j++) {
			double suma = 0;
			double sumb = 0;
			double sum2 = 0;
			for (int k : orbsOfAtom[j]) {
				suma += densityderiva.get(k, k);
				sumb += densityderivb.get(k, k);
				sum2 += densityderiv2.get(k, k);
			}

			populationsderiva[j] = -suma;
			populationsderivb[j] = -sumb;
			populationsderiv2[j] = -sum2;
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


		double[] chargedipa = new double[]{0, 0, 0};
		double[] chargedipb = new double[]{0, 0, 0};
		double[] chargedip2 = new double[]{0, 0, 0};


		for (int j = 0; j < atoms.length; j++) {
			chargedipa[0] += Constants.DIPOLECONV * populationsderiva[j] * (atoms[j].getCoordinates()[0] - com[0]);
			chargedipa[1] += Constants.DIPOLECONV * populationsderiva[j] * (atoms[j].getCoordinates()[1] - com[1]);
			chargedipa[2] += Constants.DIPOLECONV * populationsderiva[j] * (atoms[j].getCoordinates()[2] - com[2]);

			chargedipb[0] += Constants.DIPOLECONV * populationsderivb[j] * (atoms[j].getCoordinates()[0] - com[0]);
			chargedipb[1] += Constants.DIPOLECONV * populationsderivb[j] * (atoms[j].getCoordinates()[1] - com[1]);
			chargedipb[2] += Constants.DIPOLECONV * populationsderivb[j] * (atoms[j].getCoordinates()[2] - com[2]);

			chargedip2[0] += Constants.DIPOLECONV * populationsderiv2[j] * (atoms[j].getCoordinates()[0] - com[0]);
			chargedip2[1] += Constants.DIPOLECONV * populationsderiv2[j] * (atoms[j].getCoordinates()[1] - com[1]);
			chargedip2[2] += Constants.DIPOLECONV * populationsderiv2[j] * (atoms[j].getCoordinates()[2] - com[2]);
		}


		double[] hybriddipa = new double[]{0, 0, 0};
		double[] hybriddipb = new double[]{0, 0, 0};
		double[] hybriddip2 = new double[]{0, 0, 0};


		for (int j = 0; j < atoms.length; j++) {

			if (orbsOfAtom[j].length > 1) {//exclude hydrogen
				hybriddipa[0] = hybriddipa[0] - Constants.DIPOLECONV * 2 * atoms[j].D1() *
						densityderiva.get(orbsOfAtom[j][0], orbsOfAtom[j][1]);
				hybriddipa[1] = hybriddipa[1] - Constants.DIPOLECONV * 2 * atoms[j].D1() *
						densityderiva.get(orbsOfAtom[j][0], orbsOfAtom[j][2]);
				hybriddipa[2] = hybriddipa[2] - Constants.DIPOLECONV * 2 * atoms[j].D1() *
						densityderiva.get(orbsOfAtom[j][0], orbsOfAtom[j][3]);

				if (atoms[j].getAtomProperties().getZ() == Z1) {
					hybriddipa[0] = hybriddipa[0] -
							Constants.DIPOLECONV * 2 * D1deriva * densityMatrix.get(orbsOfAtom[j][0],
									orbsOfAtom[j][1]);
					hybriddipa[1] = hybriddipa[1] -
							Constants.DIPOLECONV * 2 * D1deriva * densityMatrix.get(orbsOfAtom[j][0],
									orbsOfAtom[j][2]);
					hybriddipa[2] = hybriddipa[2] -
							Constants.DIPOLECONV * 2 * D1deriva * densityMatrix.get(orbsOfAtom[j][0],
									orbsOfAtom[j][3]);
				}

				hybriddipb[0] = hybriddipb[0] - Constants.DIPOLECONV * 2 * atoms[j].D1() *
						densityderivb.get(orbsOfAtom[j][0], orbsOfAtom[j][1]);
				hybriddipb[1] = hybriddipb[1] - Constants.DIPOLECONV * 2 * atoms[j].D1() *
						densityderivb.get(orbsOfAtom[j][0], orbsOfAtom[j][2]);
				hybriddipb[2] = hybriddipb[2] - Constants.DIPOLECONV * 2 * atoms[j].D1() *
						densityderivb.get(orbsOfAtom[j][0], orbsOfAtom[j][3]);

				if (atoms[j].getAtomProperties().getZ() == Z2) {
					hybriddipb[0] = hybriddipb[0] -
							Constants.DIPOLECONV * 2 * D1derivb * densityMatrix.get(orbsOfAtom[j][0],
									orbsOfAtom[j][1]);
					hybriddipb[1] = hybriddipb[1] -
							Constants.DIPOLECONV * 2 * D1derivb * densityMatrix.get(orbsOfAtom[j][0],
									orbsOfAtom[j][2]);
					hybriddipb[2] = hybriddipb[2] -
							Constants.DIPOLECONV * 2 * D1derivb * densityMatrix.get(orbsOfAtom[j][0],
									orbsOfAtom[j][3]);
				}

				hybriddip2[0] = hybriddip2[0] - Constants.DIPOLECONV * 2 * atoms[j].D1() *
						densityderiv2.get(orbsOfAtom[j][0], orbsOfAtom[j][1]);
				hybriddip2[1] = hybriddip2[1] - Constants.DIPOLECONV * 2 * atoms[j].D1() *
						densityderiv2.get(orbsOfAtom[j][0], orbsOfAtom[j][2]);
				hybriddip2[2] = hybriddip2[2] - Constants.DIPOLECONV * 2 * atoms[j].D1() *
						densityderiv2.get(orbsOfAtom[j][0], orbsOfAtom[j][3]);

				if (atoms[j].getAtomProperties().getZ() == Z1) {
					hybriddip2[0] = hybriddip2[0] -
							Constants.DIPOLECONV * 2 * D1deriva * densityderivb.get(orbsOfAtom[j][0],
									orbsOfAtom[j][1]);
					hybriddip2[1] = hybriddip2[1] -
							Constants.DIPOLECONV * 2 * D1deriva * densityderivb.get(orbsOfAtom[j][0],
									orbsOfAtom[j][2]);
					hybriddip2[2] = hybriddip2[2] -
							Constants.DIPOLECONV * 2 * D1deriva * densityderivb.get(orbsOfAtom[j][0],
									orbsOfAtom[j][3]);
				}

				if (atoms[j].getAtomProperties().getZ() == Z2) {
					hybriddip2[0] = hybriddip2[0] -
							Constants.DIPOLECONV * 2 * D1derivb * densityderiva.get(orbsOfAtom[j][0],
									orbsOfAtom[j][1]);
					hybriddip2[1] = hybriddip2[1] -
							Constants.DIPOLECONV * 2 * D1derivb * densityderiva.get(orbsOfAtom[j][0],
									orbsOfAtom[j][2]);
					hybriddip2[2] = hybriddip2[2] -
							Constants.DIPOLECONV * 2 * D1derivb * densityderiva.get(orbsOfAtom[j][0],
									orbsOfAtom[j][3]);
				}

				if (atoms[j].getAtomProperties().getZ() == Z1 && D1deriv2 != 0) {
					hybriddip2[0] = hybriddip2[0] -
							Constants.DIPOLECONV * 2 * D1deriv2 * densityMatrix.get(orbsOfAtom[j][0],
									orbsOfAtom[j][1]);
					hybriddip2[1] = hybriddip2[1] -
							Constants.DIPOLECONV * 2 * D1deriv2 * densityMatrix.get(orbsOfAtom[j][0],
									orbsOfAtom[j][2]);
					hybriddip2[2] = hybriddip2[2] -
							Constants.DIPOLECONV * 2 * D1deriv2 * densityMatrix.get(orbsOfAtom[j][0],
									orbsOfAtom[j][3]);
				}
			}
		}

		double[] dipoletota = new double[]{chargedipa[0] + hybriddipa[0], chargedipa[1] + hybriddipa[1],
				chargedipa[2] + hybriddipa[2]};

		double[] dipoletotb = new double[]{chargedipb[0] + hybriddipb[0], chargedipb[1] + hybriddipb[1],
				chargedipb[2] + hybriddipb[2]};

		double[] dipoletot2 = new double[]{chargedip2[0] + hybriddip2[0], chargedip2[1] + hybriddip2[1],
				chargedip2[2] + hybriddip2[2]};


		return (dipoletot2[0] * soln.dipoletot[0] + dipoletota[0] * dipoletotb[0] + dipoletot2[1] * soln.dipoletot[1] +
				dipoletota[1] * dipoletotb[1] + dipoletot2[2] * soln.dipoletot[2] + dipoletota[2] * dipoletotb[2]) /
				soln.dipole - (dipoletota[0] * soln.dipoletot[0] + dipoletota[1] * soln.dipoletot[1] +
				dipoletota[2] * soln.dipoletot[2]) *
				(dipoletotb[0] * soln.dipoletot[0] + dipoletotb[1] * soln.dipoletot[1] +
						dipoletotb[2] * soln.dipoletot[2]) / Pow.pow(soln.dipole, 3);
	}

	public static double homoDeriv2(SolutionR soln, SimpleMatrix xA, SimpleMatrix xB, SimpleMatrix totalderiv,
									SimpleMatrix Fderiva, SimpleMatrix Fderivb, SimpleMatrix Fderiv2) {


		SimpleMatrix gammadiag = SimpleMatrix.diag(xA.mult(xB).plus(xB.mult(xA)).diag().getDDRM().data).scale(0.5);

		SimpleMatrix gammaremainder = gammaMatrix(soln, totalderiv);

		SimpleMatrix gamma = gammadiag.plus(gammaremainder);

		SimpleMatrix F = soln.F;

		SimpleMatrix C = soln.C;

		SimpleMatrix Cderiva = C.mult(xA);

		SimpleMatrix Cderivb = C.mult(xB);

		SimpleMatrix Cderiv2 = C.mult(gamma);

		SimpleMatrix E = Cderiv2.transpose().mult(F.mult(C)).plus(Cderiva.transpose().mult(Fderivb.mult(C)))
				.plus(Cderiva.transpose().mult(F.mult(Cderivb)));

		E = E.plus(Cderivb.transpose().mult(Fderiva.mult(C))).plus(C.transpose().mult(Fderiv2.mult(C)))
				.plus(C.transpose().mult(Fderiva.mult(Cderivb)));

		E = E.plus(Cderivb.transpose().mult(F.mult(Cderiva))).plus(C.transpose().mult(Fderivb.mult(Cderiva)))
				.plus(C.transpose().mult(F.mult(Cderiv2)));

		return E.diag().get((int) (soln.nElectrons / 2.0) - 1);
	}

	public static double homoDeriv2(SolutionU soln, SimpleMatrix xAalpha, SimpleMatrix xBalpha,
									SimpleMatrix totalderivalpha, SimpleMatrix Faderiva, SimpleMatrix Faderivb,
									SimpleMatrix Faderiv2) {


		SimpleMatrix gammadiag =
				SimpleMatrix.diag(xAalpha.mult(xBalpha).plus(xBalpha.mult(xAalpha)).diag().getDDRM().data).scale(0.5);

		SimpleMatrix gammaremainder = gammaMatrix(soln, totalderivalpha);

		SimpleMatrix gamma = gammadiag.plus(gammaremainder);

		SimpleMatrix Fa = soln.Fa;

		SimpleMatrix Ca = soln.Cta.transpose();

		SimpleMatrix Caderiva = Ca.mult(xAalpha);

		SimpleMatrix Caderivb = Ca.mult(xBalpha);

		SimpleMatrix Caderiv2 = Ca.mult(gamma);

		SimpleMatrix Ea = Caderiv2.transpose().mult(Fa.mult(Ca)).plus(Caderiva.transpose().mult(Faderivb.mult(Ca)))
				.plus(Caderiva.transpose().mult(Fa.mult(Caderivb)));

		Ea = Ea.plus(Caderivb.transpose().mult(Faderiva.mult(Ca))).plus(Ca.transpose().mult(Faderiv2.mult(Ca)))
				.plus(Ca.transpose().mult(Faderiva.mult(Caderivb)));

		Ea = Ea.plus(Caderivb.transpose().mult(Fa.mult(Caderiva))).plus(Ca.transpose().mult(Faderivb.mult(Caderiva)))
				.plus(Ca.transpose().mult(Fa.mult(Caderiv2)));

		return Ea.diag().get(soln.rm.nOccAlpha - 1);
	}


	@Deprecated
	public static SimpleMatrix staticMatrix(SolutionR soln, SimpleMatrix Fstatictotal, SimpleMatrix FstaticA,
											SimpleMatrix FstaticB, SimpleMatrix xvectorA, SimpleMatrix xvectorB,
											int Z1,
											int param1, int Z2, int param2) {

		SimpleMatrix densityA = PopleThiel.densityDeriv(soln, xvectorA);
		SimpleMatrix densityB = PopleThiel.densityDeriv(soln, xvectorB);

		SimpleMatrix GaB = ParamSecondDerivative.Gderivstatic(soln, densityB, Z1, param1);
		SimpleMatrix GbA = ParamSecondDerivative.Gderivstatic(soln, densityA, Z2, param2);

		SimpleMatrix FockA = FstaticA.plus(PopleThiel.responseMatrix(soln, densityA));
		SimpleMatrix FockB = FstaticB.plus(PopleThiel.responseMatrix(soln, densityB));

		SimpleMatrix FA = soln.Ct.mult(FockA.mult(soln.C));
		SimpleMatrix diagFA = SimpleMatrix.diag(FA.diag().getDDRM().data);
		SimpleMatrix FB = soln.Ct.mult(FockB.mult(soln.C));
		SimpleMatrix diagFB = SimpleMatrix.diag(FB.diag().getDDRM().data);


		SimpleMatrix xA = ParamDerivative.xMatrix(soln, FA);
		SimpleMatrix xB = ParamDerivative.xMatrix(soln, FB);

		SimpleMatrix omega = PopleThiel.responseMatrix(soln, ParamSecondDerivative.densityDeriv2static(soln, xA, xB));

		SimpleMatrix Phi = Fstatictotal.plus(GaB).plus(GbA).plus(omega);

		SimpleMatrix matrix = soln.Ct.mult(Phi).mult(soln.C);

		matrix = matrix.plus(xB.transpose().mult(diagFA)).plus(diagFA.mult(xB)).plus(xA.transpose().mult(diagFB))
				.plus(diagFB.mult(xA));
		matrix = matrix.minus(xA.transpose().mult(SimpleMatrix.diag(soln.E.getDDRM().data)).mult(xB))
				.minus(xB.transpose().mult(SimpleMatrix.diag(soln.E.getDDRM().data)).mult(xA));
		matrix = matrix.minus(SimpleMatrix.diag(soln.E.getDDRM().data).mult(xA).mult(xB))
				.minus(SimpleMatrix.diag(soln.E.getDDRM().data).mult(xB).mult(xA));

		return matrix;
	}

	@Deprecated
	public static SimpleMatrix[] staticMatrix(SolutionU soln, SimpleMatrix Fstatictotalalpha,
											  SimpleMatrix Fstatictotalbeta, SimpleMatrix FstaticAalpha,
											  SimpleMatrix FstaticAbeta, SimpleMatrix FstaticBalpha,
											  SimpleMatrix FstaticBbeta, SimpleMatrix xvectorA, SimpleMatrix xvectorB,
											  int Z1, int param1, int Z2, int param2) {

		SimpleMatrix[] densitiesA = PopleThiel.densityDeriv(soln, xvectorA);

		SimpleMatrix densityAalpha = densitiesA[0];
		SimpleMatrix densityAbeta = densitiesA[1];

		SimpleMatrix[] densitiesB = PopleThiel.densityDeriv(soln, xvectorB);

		SimpleMatrix densityBalpha = densitiesB[0];
		SimpleMatrix densityBbeta = densitiesB[1];

		SimpleMatrix[] GaB = ParamSecondDerivative.Gderivstatic(soln, densityBalpha, densityBbeta, Z1, param1);

		SimpleMatrix GaBalpha = GaB[0];

		SimpleMatrix GaBbeta = GaB[1];

		SimpleMatrix[] GbA = ParamSecondDerivative.Gderivstatic(soln, densityAalpha, densityAbeta, Z2, param2);

		SimpleMatrix GbAalpha = GbA[0];

		SimpleMatrix GbAbeta = GbA[1];

		SimpleMatrix[] Ra = PopleThiel.responseMatrix(soln, densitiesA);

		SimpleMatrix FockAalpha = FstaticAalpha.plus(Ra[0]);
		SimpleMatrix FockAbeta = FstaticAbeta.plus(Ra[1]);

		SimpleMatrix[] Rb = PopleThiel.responseMatrix(soln, densitiesB);

		SimpleMatrix FockBalpha = FstaticBalpha.plus(Rb[0]);
		SimpleMatrix FockBbeta = FstaticBbeta.plus(Rb[1]);

		SimpleMatrix FAalpha = soln.Cta.mult(FockAalpha).mult(soln.Cta.transpose());
		SimpleMatrix FAbeta = soln.Ctb.mult(FockAbeta).mult(soln.Ctb.transpose());

		SimpleMatrix FBalpha = soln.Cta.mult(FockBalpha).mult(soln.Cta.transpose());
		SimpleMatrix FBbeta = soln.Ctb.mult(FockBbeta).mult(soln.Ctb.transpose());

		SimpleMatrix diagFAalpha = SimpleMatrix.diag(FAalpha.diag().getDDRM().data);
		SimpleMatrix diagFBalpha = SimpleMatrix.diag(FBalpha.diag().getDDRM().data);

		SimpleMatrix diagFAbeta = SimpleMatrix.diag(FAbeta.diag().getDDRM().data);
		SimpleMatrix diagFBbeta = SimpleMatrix.diag(FBbeta.diag().getDDRM().data);

		SimpleMatrix[] xA = ParamDerivative.xMatrix(soln, FAalpha, FAbeta);

		SimpleMatrix xAalpha = xA[0];
		SimpleMatrix xAbeta = xA[1];

		SimpleMatrix[] xB = ParamDerivative.xMatrix(soln, FBalpha, FBbeta);

		SimpleMatrix xBalpha = xB[0];
		SimpleMatrix xBbeta = xB[1];

		SimpleMatrix[] omega = PopleThiel.responseMatrix(soln, ParamSecondDerivative.densityDeriv2static(soln, xA,
				xB));

		SimpleMatrix omegaalpha = omega[0];

		SimpleMatrix omegabeta = omega[1];

		SimpleMatrix Phialpha = Fstatictotalalpha.plus(GaBalpha).plus(GbAalpha).plus(omegaalpha);
		SimpleMatrix Phibeta = Fstatictotalbeta.plus(GaBbeta).plus(GbAbeta).plus(omegabeta);

		SimpleMatrix matrixalpha = soln.Cta.mult(Phialpha).mult(soln.Cta.transpose());
		SimpleMatrix matrixbeta = soln.Ctb.mult(Phibeta).mult(soln.Ctb.transpose());

		matrixalpha = matrixalpha.plus(xBalpha.transpose().mult(diagFAalpha)).plus(diagFAalpha.mult(xBalpha))
				.plus(xAalpha.transpose().mult(diagFBalpha)).plus(diagFBalpha.mult(xAalpha));

		matrixalpha =
				matrixalpha.minus(xAalpha.transpose().mult(SimpleMatrix.diag(soln.Ea.getDDRM().data)).mult(xBalpha))
						.minus(xBalpha.transpose().mult(SimpleMatrix.diag(soln.Ea.getDDRM().data)).mult(xAalpha));

		matrixalpha = matrixalpha.minus(SimpleMatrix.diag(soln.Ea.getDDRM().data).mult(xAalpha).mult(xBalpha))
				.minus(SimpleMatrix.diag(soln.Ea.getDDRM().data).mult(xBalpha).mult(xAalpha));


		matrixbeta = matrixbeta.plus(xBbeta.transpose().mult(diagFAbeta)).plus(diagFAbeta.mult(xBbeta))
				.plus(xAbeta.transpose().mult(diagFBbeta)).plus(diagFBbeta.mult(xAbeta));

		matrixbeta = matrixbeta.minus(xAbeta.transpose().mult(SimpleMatrix.diag(soln.Eb.getDDRM().data)).mult(xBbeta))
				.minus(xBbeta.transpose().mult(SimpleMatrix.diag(soln.Eb.getDDRM().data)).mult(xAbeta));

		matrixbeta = matrixbeta.minus(SimpleMatrix.diag(soln.Eb.getDDRM().data).mult(xAbeta).mult(xBbeta))
				.minus(SimpleMatrix.diag(soln.Eb.getDDRM().data).mult(xBbeta).mult(xAbeta));

		return new SimpleMatrix[]{matrixalpha, matrixbeta};
	}

	@Deprecated
	public static SimpleMatrix[] staticFockDeriv(SolutionU soln, SimpleMatrix Fstatictotalalpha,
												 SimpleMatrix Fstatictotalbeta, SimpleMatrix FstaticAalpha,
												 SimpleMatrix FstaticAbeta, SimpleMatrix FstaticBalpha,
												 SimpleMatrix FstaticBbeta, SimpleMatrix xvectorA,
												 SimpleMatrix xvectorB, int Z1, int param1, int Z2, int param2) {

		SimpleMatrix[] densitiesA = PopleThiel.densityDeriv(soln, xvectorA);

		SimpleMatrix densityAalpha = densitiesA[0];
		SimpleMatrix densityAbeta = densitiesA[1];

		SimpleMatrix[] densitiesB = PopleThiel.densityDeriv(soln, xvectorB);

		SimpleMatrix densityBalpha = densitiesB[0];
		SimpleMatrix densityBbeta = densitiesB[1];

		SimpleMatrix[] GaB = ParamSecondDerivative.Gderivstatic(soln, densityBalpha, densityBbeta, Z1, param1);

		SimpleMatrix GaBalpha = GaB[0];

		SimpleMatrix GaBbeta = GaB[1];

		SimpleMatrix[] GbA = ParamSecondDerivative.Gderivstatic(soln, densityAalpha, densityAbeta, Z2, param2);

		SimpleMatrix GbAalpha = GbA[0];

		SimpleMatrix GbAbeta = GbA[1];

		SimpleMatrix[] Ra = PopleThiel.responseMatrix(soln, densitiesA);

		SimpleMatrix FockAalpha = FstaticAalpha.plus(Ra[0]);
		SimpleMatrix FockAbeta = FstaticAbeta.plus(Ra[1]);

		SimpleMatrix[] Rb = PopleThiel.responseMatrix(soln, densitiesB);

		SimpleMatrix FockBalpha = FstaticBalpha.plus(Rb[0]);
		SimpleMatrix FockBbeta = FstaticBbeta.plus(Rb[1]);

		SimpleMatrix FAalpha = soln.Cta.mult(FockAalpha).mult(soln.Cta.transpose());
		SimpleMatrix FAbeta = soln.Ctb.mult(FockAbeta).mult(soln.Ctb.transpose());

		SimpleMatrix FBalpha = soln.Cta.mult(FockBalpha).mult(soln.Cta.transpose());
		SimpleMatrix FBbeta = soln.Ctb.mult(FockBbeta).mult(soln.Ctb.transpose());

		SimpleMatrix[] xA = ParamDerivative.xMatrix(soln, FAalpha, FAbeta);

		SimpleMatrix[] xB = ParamDerivative.xMatrix(soln, FBalpha, FBbeta);

		SimpleMatrix[] omega = PopleThiel.responseMatrix(soln, ParamSecondDerivative.densityDeriv2static(soln, xA,
				xB));

		SimpleMatrix omegaalpha = omega[0];
		SimpleMatrix omegabeta = omega[1];

		SimpleMatrix Phialpha = Fstatictotalalpha.plus(GaBalpha).plus(GbAalpha).plus(omegaalpha);
		SimpleMatrix Phibeta = Fstatictotalbeta.plus(GaBbeta).plus(GbAbeta).plus(omegabeta);

		return new SimpleMatrix[]{Phialpha, Phibeta};
	}

	@Deprecated
	public static SimpleMatrix utilitylazy(SolutionR soln, SimpleMatrix Fstatictotal, SimpleMatrix FstaticA,
										   SimpleMatrix FstaticB, SimpleMatrix xvectorA, SimpleMatrix xvectorB, int Z1,
										   int param1, int Z2, int param2) {

		SimpleMatrix densityA = PopleThiel.densityDeriv(soln, xvectorA);
		SimpleMatrix densityB = PopleThiel.densityDeriv(soln, xvectorB);

		SimpleMatrix FockA = FstaticA.plus(PopleThiel.responseMatrix(soln, densityA));
		SimpleMatrix FockB = FstaticB.plus(PopleThiel.responseMatrix(soln, densityB));

		SimpleMatrix FA = soln.Ct.mult(FockA.mult(soln.C));
		SimpleMatrix FB = soln.Ct.mult(FockB.mult(soln.C));

		SimpleMatrix xA = ParamDerivative.xMatrix(soln, FA);
		SimpleMatrix xB = ParamDerivative.xMatrix(soln, FB);

		SimpleMatrix densityderivstatic = ParamSecondDerivative.densityDeriv2static(soln, xA, xB);

		return densityderivstatic;
	}

	@Deprecated
	public static SimpleMatrix[] utilitylazy(SolutionU soln, SimpleMatrix Fstatictotalalpha,
											 SimpleMatrix Fstatictotalbeta, SimpleMatrix FstaticAalpha,
											 SimpleMatrix FstaticAbeta, SimpleMatrix FstaticBalpha,
											 SimpleMatrix FstaticBbeta, SimpleMatrix xvectorA, SimpleMatrix xvectorB,
											 int Z1, int param1, int Z2, int param2) {

		SimpleMatrix[] densitiesA = PopleThiel.densityDeriv(soln, xvectorA);

		SimpleMatrix[] densitiesB = PopleThiel.densityDeriv(soln, xvectorB);

		SimpleMatrix[] Ra = PopleThiel.responseMatrix(soln, densitiesA);

		SimpleMatrix FockAalpha = FstaticAalpha.plus(Ra[0]);
		SimpleMatrix FockAbeta = FstaticAbeta.plus(Ra[1]);

		SimpleMatrix[] Rb = PopleThiel.responseMatrix(soln, densitiesB);

		SimpleMatrix FockBalpha = FstaticBalpha.plus(Rb[0]);
		SimpleMatrix FockBbeta = FstaticBbeta.plus(Rb[1]);

		SimpleMatrix FAalpha = soln.Cta.mult(FockAalpha).mult(soln.Cta.transpose());
		SimpleMatrix FAbeta = soln.Ctb.mult(FockAbeta).mult(soln.Ctb.transpose());

		SimpleMatrix FBalpha = soln.Cta.mult(FockBalpha).mult(soln.Cta.transpose());
		SimpleMatrix FBbeta = soln.Ctb.mult(FockBbeta).mult(soln.Ctb.transpose());


		SimpleMatrix[] xA = ParamDerivative.xMatrix(soln, FAalpha, FAbeta);

		SimpleMatrix[] xB = ParamDerivative.xMatrix(soln, FBalpha, FBbeta);

		return ParamSecondDerivative.densityDeriv2static(soln, xA, xB);
	}

	@Deprecated
	public static SimpleMatrix[] gammaArrayThiel(SolutionU soln, SimpleMatrix[] matrixderivalpha,
												 SimpleMatrix[] matrixderivbeta) {

		StopWatch sw = new StopWatch();
		sw.start();


		int NOccAlpha = soln.rm.nOccAlpha;
		int NOccBeta = soln.rm.nOccBeta;

		int NVirtAlpha = soln.rm.nVirtAlpha;
		int NVirtBeta = soln.rm.nVirtBeta;

		SimpleMatrix[] xarray = new SimpleMatrix[matrixderivalpha.length];

		SimpleMatrix[] rarray = new SimpleMatrix[matrixderivalpha.length];

		SimpleMatrix[] dirs = new SimpleMatrix[matrixderivalpha.length];

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

					F.set(count1, 0, matrixderivalpha[a].get(i, j + NOccAlpha));

					count1++;
				}
			}

			for (int i = 0; i < NOccBeta; i++) { // kappa
				for (int j = 0; j < NVirtBeta; j++) { // i

					F.set(count1, 0, matrixderivbeta[a].get(i, j + NOccBeta));


					count1++;
				}
			}

			F = D.mult(F);

			xarray[a] = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);

			rarray[a] = F.copy();

			dirs[a] = F.copy();
		}


		if (dirs[0].numRows() == 0) {
			SimpleMatrix[] densityderivs = new SimpleMatrix[matrixderivalpha.length];

			for (int i = 0; i < densityderivs.length; i++) {
				densityderivs[i] = new SimpleMatrix(0, 0);
			}

			return densityderivs;
		}


		while (Utils.numNotNull(rarray) > 0) {

			ArrayList<SimpleMatrix> d = new ArrayList<>();

			ArrayList<SimpleMatrix> p = new ArrayList<>();

//			System.err.println(
//					"It's still running, don't worry: " + numNotNull(rarray));

			for (int i = 0; i < rarray.length; i++) {

				if (rarray[i] != null) {

					d.add(dirs[i].copy());
					p.add(D.mult(computeResponseVectorsThiel(dirs[i], soln)));
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

					double val = p.get(j).transpose().mult(d.get(i)).get(0, 0) +
							p.get(i).transpose().mult(d.get(j)).get(0, 0);
					solver.set(i, j, val);
					solver.set(j, i, val);
				}
			}

			SimpleMatrix alpha = solver.solve(rhsvec);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {


					for (int i = 0; i < alpha.numRows(); i++) {
						xarray[a] = xarray[a].plus(d.get(i).scale(alpha.get(i, a)));
						rarray[a] = rarray[a].minus(p.get(i).scale(alpha.get(i, a)));

					}

					if (Utils.mag(rarray[a]) < 1E-8) {//todo change this
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
						rhs.set(i, 0, -rarray[a].transpose().mult(p.get(i)).get(0, 0));

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

	public static boolean verifyEquations(SolutionR soln, int Z1, int param1, int Z2, int param2) {

		SolutionR solnprime = (SolutionR) soln.withNewAtoms(Utils.perturbAtomParams(soln.atoms, Z1, param1));

		solnprime.compute();

		SimpleMatrix H1 = ParamDerivative.staticDeriv(soln, Z1, 0)[0][param1];
		SimpleMatrix H2 = ParamDerivative.staticDeriv(soln, Z2, 0)[0][param2];


		SimpleMatrix F2 = ParamDerivative.staticDeriv(soln, Z2, 0)[1][param2];
		SimpleMatrix F2prime = ParamDerivative.staticDeriv(solnprime, Z2, 0)[1][param2];
		SimpleMatrix F1 = ParamDerivative.staticDeriv(soln, Z1, 0)[1][param1];

		SimpleMatrix x2 = PopleThiel.pt(soln, PopleThiel.toMO(soln.CtOcc, soln.CVirt, new SimpleMatrix[]{F2}))[0];
		SimpleMatrix x2prime = PopleThiel.pt(solnprime,
				PopleThiel.toMO(solnprime.CtOcc, solnprime.CVirt, new SimpleMatrix[]{F2prime}))[0];
		SimpleMatrix x1 = PopleThiel.pt(soln, PopleThiel.toMO(soln.CtOcc, soln.CVirt, new SimpleMatrix[]{F1}))[0];

		SimpleMatrix D1 = PopleThiel.densityDeriv(soln, x1);

		SimpleMatrix R1 = PopleThiel.responseMatrix(soln, D1);

		SimpleMatrix C1 = soln.C.mult(ParamDerivative.xMatrix(soln, soln.Ct.mult(F1.plus(R1)).mult(soln.C)));


		SimpleMatrix D2prime = PopleThiel.densityDeriv(solnprime, x2prime);


		SimpleMatrix D2 = PopleThiel.densityDeriv(soln, x2);

		SimpleMatrix R2 = PopleThiel.responseMatrix(soln, D2);

		SimpleMatrix C2 = soln.C.mult(ParamDerivative.xMatrix(soln, soln.Ct.mult(F2.plus(R2)).mult(soln.C)));

		SimpleMatrix densityderiv2finite = D2prime.minus(D2).scale(1 / Constants.LAMBDA);


		// second order Fstatic
		SimpleMatrix Fstatictotal = ParamSecondDerivative.Hderiv2(soln, Z1, param1, Z2, param2)
				.plus(ParamSecondDerivative.Gderiv2static(soln, Z1, param1, Z2, param2));

		SimpleMatrix Hstatic = ParamSecondDerivative.Hderiv2(soln, Z1, param1, Z2, param2);


		SimpleMatrix densityderiv2static =
				ParamSecondDerivative.utilitylazy(soln, Fstatictotal, F1, F2, x1, x2, Z1, param1, Z2, param2);

		SimpleMatrix rhsmat =
				ParamSecondDerivative.staticMatrix(soln, Fstatictotal, F1, F2, x1, x2, Z1, param1, Z2, param2);

//		SimpleMatrix gammavec = PopleThiel.pople(soln, new SimpleMatrix[]{rhsmat})[0];

		SimpleMatrix rhsmat2 = rhsmat.extractMatrix(0, soln.rm.nOccAlpha, soln.rm.nOccAlpha, soln.rm.nOrbitals);

		SimpleMatrix newgammavec = PopleThiel.pt(soln, new SimpleMatrix[]{rhsmat2})[0];

		SimpleMatrix densityderiv2response = PopleThiel.densityDeriv(soln, newgammavec);

		SimpleMatrix densityderiv2 = densityderiv2response.plus(densityderiv2static);

		double Hfderiv = 1 / Constants.LAMBDA *
				(ParamDerivative.HfDeriv(solnprime, Z2, param2) - ParamDerivative.HfDeriv(soln, Z2, param2));

		double Hfderivtest = ParamSecondDerivative.HfDeriv2(soln, Z1, param1, Z2, param2, H1, F1, D2, 0);

		double Hfderivtest2 = ParamSecondDerivative.HfDeriv2(soln, Z1, param1, Z2, param2, H2, F2, D1, 1);

		SimpleMatrix x1mat =
				ParamDerivative.xMatrix(soln, soln.Ct.mult(F1.plus(PopleThiel.responseMatrix(soln, D1))).mult(soln.C));

		SimpleMatrix x2mat =
				ParamDerivative.xMatrix(soln, soln.Ct.mult(F2.plus(PopleThiel.responseMatrix(soln, D2))).mult(soln.C));

		SimpleMatrix x2matprime = ParamDerivative.xMatrix(solnprime,
				solnprime.Ct.mult(F2prime.plus(PopleThiel.responseMatrix(solnprime, D2prime))).mult(solnprime.C));

		SimpleMatrix totalderiv =
				rhsmat.plus(soln.Ct.mult(PopleThiel.responseMatrix(soln, densityderiv2response)).mult(soln.C));

		SimpleMatrix Fderiv2 =
				staticFockDeriv(soln, Fstatictotal, D1, D2, densityderiv2static, Z1, param1, Z2, param2).plus(
						PopleThiel.responseMatrix(soln, densityderiv2response));

		SimpleMatrix Fderiva = F1.plus(PopleThiel.responseMatrix(soln, D1));

		SimpleMatrix Fderivb = F2.plus(PopleThiel.responseMatrix(soln, D2));

		SimpleMatrix Fderivbprime = F2prime.plus(PopleThiel.responseMatrix(solnprime, D2prime));

		double IEtest = homoDeriv2(soln, x1mat, x2mat, totalderiv, Fderiva, Fderivb, Fderiv2);

		double IEderiv = ParamDerivative.homoDeriv(soln, x2mat, soln.Ct.mult(Fderivb).mult(soln.C));

		double IEderivprime =
				ParamDerivative.homoDeriv(solnprime, x2matprime, solnprime.Ct.mult(Fderivbprime).mult(solnprime.C));


		double IEderiv2 = 1 / Constants.LAMBDA * (IEderivprime - IEderiv);

		System.err.println("----TESTING IE DERIVATIVES----");


		System.err.println(IEderiv2);
		System.err.println(IEtest);

		System.err.println("----TESTING HF DERIVATIVES----");


		//System.err.println (densityderiv2);

		//System.err.println (densityderiv2finite);

		System.err.println(Hfderiv);
		System.err.println(Hfderivtest);
		System.err.println(Hfderivtest2);

		double dipolederivtest = dipoleDeriv2(soln, D1, D2, densityderiv2, Z1, param1, Z2, param2);

		double dipolederiv = ParamDerivative.nddoDipoleDeriv(soln, D2, Z2, param2);

		double dipolederivprime = ParamDerivative.nddoDipoleDeriv(solnprime, D2prime, Z2, param2);

		double dipolederiv2 = 1 / Constants.LAMBDA * (dipolederivprime - dipolederiv);

		System.err.println("----TESTING DIPOLE DERIVATIVES----");


		System.err.println(dipolederiv2);
		System.err.println(dipolederivtest);


		System.err.println("----TESTING EXPGEOM DERIVATIVES----");

		double gradderiv = ParamGeometryDerivative.gradDeriv(soln, 0, 0, Z2, param2, D2);

		double gradderivprime = ParamGeometryDerivative.gradDeriv(solnprime, 0, 0, Z2, param2, D2prime);

		double gradderiv2test =
				ParamGeometrySecondDerivative.gradderiv2(soln, 0, 0, Z1, param1, Z2, param2, D1, D2, densityderiv2);
		double gradderiv2 = 1 / Constants.LAMBDA * (gradderivprime - gradderiv);

		System.err.println(gradderiv2);


		System.err.println(gradderiv2test);

		System.err.println("----TESTING MY CLOWNERY----");

		SimpleMatrix geomGradVector = new SimpleMatrix(soln.atoms.length * 3, 1);

		for (int i = 0, count = 0; i < soln.atoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				double d = GeometryDerivative.grad(soln, i, j);
				geomGradVector.set(count, d);
				count++;
			}
		}

		SimpleMatrix geomGradVectorprime = new SimpleMatrix(solnprime.atoms.length * 3, 1);

		for (int i = 0, count = 0; i < solnprime.atoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				double d = GeometryDerivative.grad(solnprime, i, j);
				geomGradVectorprime.set(count, d);
				count++;
			}
		}

		double deriv = (geomGradVectorprime.normF() - geomGradVector.normF()) / Constants.LAMBDA;

		SimpleMatrix geomGradVectorDeriv = new SimpleMatrix(soln.atoms.length * 3, 1);

		for (int i = 0, count = 0; i < soln.atoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				double d = ParamGeometryDerivative.gradDeriv(soln, i, j, Z1, param1, D1);
				geomGradVectorDeriv.set(count, d);
				count++;
			}
		}

		double derivtest = geomGradVectorDeriv.dot(geomGradVector) / geomGradVector.normF();


		SimpleMatrix geomGradVectorDeriv2 = new SimpleMatrix(soln.atoms.length * 3, 1);

		for (int i = 0, count = 0; i < soln.atoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				double d = ParamGeometryDerivative.gradDeriv(soln, i, j, Z2, param2, D2);
				geomGradVectorDeriv2.set(count, d);
				count++;
			}
		}

		SimpleMatrix geomGradVectorDerivprime2 = new SimpleMatrix(soln.atoms.length * 3, 1);

		for (int i = 0, count = 0; i < soln.atoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				double d = ParamGeometryDerivative.gradDeriv(solnprime, i, j, Z2, param2, D2);
				geomGradVectorDerivprime2.set(count, d);
				count++;
			}
		}

		SimpleMatrix geomGradVectorSecondDeriv = new SimpleMatrix(soln.atoms.length * 3, 1);

		for (int i = 0, count = 0; i < soln.atoms.length; i++) {
			for (int j = 0; j < 3; j++) {
				double d = ParamGeometrySecondDerivative.gradderiv2(soln, i, j, Z1, param1, Z2, param2, D1, D2,
						densityderiv2);
				geomGradVectorSecondDeriv.set(count, d);
				count++;
			}
		}

		double derivtestb = geomGradVectorDeriv2.dot(geomGradVector) / geomGradVector.normF();


		double deriv2 = (geomGradVectorDerivprime2.normF() - geomGradVectorDeriv2.normF()) / Constants.LAMBDA;

		double deriv2test =
				(geomGradVectorSecondDeriv.dot(geomGradVector) + geomGradVectorDeriv.dot(geomGradVectorDeriv2) -
						derivtest * derivtestb) / geomGradVector.normF();


		System.err.println(deriv + "; " + derivtest);

		System.err.println(deriv2 + "; " + deriv2test);


		System.err.println("----TESTING CONCLUDED----");

		return Utils.isSimilar(densityderiv2, densityderiv2finite, 1E-5) && Math.abs(Hfderiv - Hfderivtest) < 1E-3 &&
				Math.abs(IEderiv2 - IEtest) < 1E-3 && Math.abs(dipolederiv2 - dipolederivtest) < 1E-3 &&
				Math.abs(gradderiv2test - gradderiv2) < 1E-5;
	}

	public static boolean verifyEquations(SolutionU soln, int Z1, int param1, int Z2, int param2) {

		SolutionU solnprime = (SolutionU) soln.withNewAtoms(Utils.perturbAtomParams(soln.atoms, Z1, param1));

		solnprime.compute();

		SimpleMatrix H1 = ParamDerivative.staticDeriv(soln, Z1, 0)[0][param1];
		SimpleMatrix H2 = ParamDerivative.staticDeriv(soln, Z2, 0)[0][param2];


		SimpleMatrix F2alpha = ParamDerivative.staticDeriv(soln, Z2, 0)[1][param2];
		SimpleMatrix F2beta = ParamDerivative.staticDeriv(soln, Z2, 0)[2][param2];

		SimpleMatrix F2alphaprime = ParamDerivative.staticDeriv(solnprime, Z2, 0)[1][param2];
		SimpleMatrix F2betaprime = ParamDerivative.staticDeriv(solnprime, Z2, 0)[2][param2];

		SimpleMatrix F1alpha = ParamDerivative.staticDeriv(soln, Z1, 0)[1][param1];
		SimpleMatrix F1beta = ParamDerivative.staticDeriv(soln, Z1, 0)[2][param1];


		SimpleMatrix x2 = PopleThiel.pt(soln, PopleThiel.toMO(soln.CtaOcc, soln.CaVirt, new SimpleMatrix[]{F2alpha}),
				PopleThiel.toMO(soln.CtbOcc, soln.CbVirt, new SimpleMatrix[]{F2beta}))[0];
		SimpleMatrix x2prime = PopleThiel.pt(solnprime,
				PopleThiel.toMO(solnprime.CtaOcc, solnprime.CaVirt, new SimpleMatrix[]{F2alphaprime}),
				PopleThiel.toMO(solnprime.CtbOcc, solnprime.CbVirt, new SimpleMatrix[]{F2betaprime}))[0];
		SimpleMatrix x1 = PopleThiel.pt(soln, PopleThiel.toMO(soln.CtaOcc, soln.CaVirt, new SimpleMatrix[]{F1alpha}),
				PopleThiel.toMO(soln.CtbOcc, soln.CbVirt, new SimpleMatrix[]{F1beta}))[0];

		SimpleMatrix[] D1 = PopleThiel.densityDeriv(soln, x1);

		SimpleMatrix[] R1 = PopleThiel.responseMatrix(soln, D1);

		SimpleMatrix[] x1matrix =
				ParamDerivative.xMatrix(soln, soln.Cta.mult(F1alpha.plus(R1[0])).mult(soln.Cta.transpose()),
						soln.Ctb.mult(F1beta.plus(R1[1])).mult(soln.Ctb.transpose()));

		SimpleMatrix C1alpha = soln.Cta.transpose().mult(x1matrix[0]);
		SimpleMatrix C1beta = soln.Ctb.transpose().mult(x1matrix[1]);

		SimpleMatrix[] D2prime = PopleThiel.densityDeriv(solnprime, x2prime);

		SimpleMatrix[] D2 = PopleThiel.densityDeriv(soln, x2);

		SimpleMatrix[] R2 = PopleThiel.responseMatrix(soln, D2);

		SimpleMatrix[] R2prime = PopleThiel.responseMatrix(solnprime, D2prime);


		SimpleMatrix[] x2matrix =
				ParamDerivative.xMatrix(soln, soln.Cta.mult(F2alpha.plus(R2[0])).mult(soln.Cta.transpose()),
						soln.Ctb.mult(F2beta.plus(R2[1])).mult(soln.Ctb.transpose()));

		SimpleMatrix[] x2primematrix = ParamDerivative.xMatrix(solnprime,
				solnprime.Cta.mult(F2alphaprime.plus(R2prime[0])).mult(solnprime.Cta.transpose()),
				solnprime.Ctb.mult(F2betaprime.plus(R2prime[1])).mult(solnprime.Ctb.transpose()));


		SimpleMatrix C2alpha = soln.Cta.transpose().mult(x2matrix[0]);
		SimpleMatrix C2beta = soln.Ctb.transpose().mult(x2matrix[1]);

		SimpleMatrix densityderiv2alphafinite = D2prime[0].minus(D2[0]).scale(1 / Constants.LAMBDA);
		SimpleMatrix densityderiv2betafinite = D2prime[1].minus(D2[1]).scale(1 / Constants.LAMBDA);

		SimpleMatrix[] Gmatrices = ParamSecondDerivative.Gderiv2static(soln, Z1, param1, Z2, param2);


		SimpleMatrix Fstatictotalalpha =
				ParamSecondDerivative.Hderiv2(soln, Z1, param1, Z2, param2).plus(Gmatrices[0]);

		SimpleMatrix Fstatictotalbeta = ParamSecondDerivative.Hderiv2(soln, Z1, param1, Z2, param2).plus(Gmatrices[1]);


		SimpleMatrix[] densityderiv2static =
				ParamSecondDerivative.utilitylazy(soln, Fstatictotalalpha, Fstatictotalbeta, F1alpha, F1beta, F2alpha,
						F2beta, x1, x2, Z1, param1, Z2, param2);

		SimpleMatrix[] rhsmat =
				ParamSecondDerivative.staticMatrix(soln, Fstatictotalalpha, Fstatictotalbeta, F1alpha, F1beta, F2alpha,
						F2beta, x1, x2, Z1, param1, Z2, param2);

//		System.out.println("rhsmat = " + rhsmat[0]);

		SimpleMatrix gammavec = gammaArrayThiel(soln, new SimpleMatrix[]{rhsmat[0]}, new SimpleMatrix[]{rhsmat[1]})[0];

		SimpleMatrix[] densityderiv2response = PopleThiel.densityDeriv(soln, gammavec);

		SimpleMatrix densityderiv2alpha = densityderiv2response[0].plus(densityderiv2static[0]);

		SimpleMatrix densityderiv2beta = densityderiv2response[1].plus(densityderiv2static[1]);

//		System.err.println(densityderiv2alpha);

//		System.err.println(densityderiv2alphafinite);

		System.err.println("---");

//		System.err.println(densityderiv2beta);

//		System.err.println(densityderiv2betafinite);

		double Hfderiv = 1 / Constants.LAMBDA *
				(ParamDerivative.HfDeriv(solnprime, Z2, param2) - ParamDerivative.HfDeriv(soln, Z2, param2));

		double Hfderivtest =
				ParamSecondDerivative.HfDeriv2(soln, Z1, param1, Z2, param2, H1, F1alpha, F1beta, D2[0], D2[1], 0);

		double Hfderivtest2 =
				ParamSecondDerivative.HfDeriv2(soln, Z1, param1, Z2, param2, H2, F2alpha, F2beta, D1[0], D1[1], 1);

		SimpleMatrix totalderiv = rhsmat[0].plus(
				soln.Cta.mult(PopleThiel.responseMatrix(soln, densityderiv2response)[0]).mult(soln.Cta.transpose()));

		SimpleMatrix Fderiv2 =
				staticFockDeriv(soln, Fstatictotalalpha, Fstatictotalbeta, F1alpha, F1beta, F2alpha, F2beta, x1, x2,
						Z1,
						param1, Z2, param2)[0].plus(PopleThiel.responseMatrix(soln, densityderiv2response)[0]);

		SimpleMatrix Fderiva = F1alpha.plus(R1[0]);

		SimpleMatrix Fderivb = F2alpha.plus(R2[0]);

		SimpleMatrix Fderivbprime = F2alphaprime.plus(PopleThiel.responseMatrix(solnprime, D2prime)[0]);
//
		double IEtest = homoDeriv2(soln, x1matrix[0], x2matrix[0], totalderiv, Fderiva, Fderivb, Fderiv2);


		double IEderiv = ParamDerivative.homoDeriv(soln, x2matrix[0], soln.Cta.mult(Fderivb).mult(soln.Ca));
		double IEderivprime = ParamDerivative.homoDeriv(solnprime, x2primematrix[0],
				solnprime.Cta.mult(Fderivbprime).mult(solnprime.Ca));
		double IEderiv2 = 1 / Constants.LAMBDA * (IEderivprime - IEderiv);
//
		System.err.println("----TESTING IE DERIVATIVES----");
//
//
		System.err.println(IEderiv2);
		System.err.println(IEtest);
//
		System.err.println("----TESTING HF DERIVATIVES----");
//

//
		System.err.println(Hfderiv);
		System.err.println(Hfderivtest);
		System.err.println(Hfderivtest2);
//
		double dipolederivtest =
				dipoleDeriv2(soln, D1[0].plus(D1[1]), D2[0].plus(D2[1]), densityderiv2alpha.plus(densityderiv2beta),
						Z1,
						param1, Z2, param2);
//
		double dipolederiv = ParamDerivative.nddoDipoleDeriv(soln, D2[0].plus(D2[1]), Z2, param2);
//
		double dipolederivprime = ParamDerivative.nddoDipoleDeriv(solnprime, D2prime[0].plus(D2prime[1]), Z2, param2);
//
		double dipolederiv2 = 1 / Constants.LAMBDA * (dipolederivprime - dipolederiv);
//
		System.err.println("----TESTING DIPOLE DERIVATIVES----");
//
//
		System.err.println(dipolederiv2);
		System.err.println(dipolederivtest);

//
//
		System.err.println("----TESTING EXPGEOM DERIVATIVES----");
//
		double gradderiv = ParamGeometryDerivative.gradDeriv(soln, 0, 0, Z2, param2, D2[0], D2[1]);
//
		double gradderivprime = ParamGeometryDerivative.gradDeriv(solnprime, 0, 0, Z2, param2, D2prime[0], D2prime[1]);
//
		double gradderiv2 =
				ParamGeometrySecondDerivative.gradderiv2(soln, 0, 0, Z1, param1, Z2, param2, D1[0], D1[1], D2[0],
						D2[1],
						densityderiv2alpha, densityderiv2beta);
//
		System.err.println("analyt " + gradderiv2);
//
		double gradderiv2test = 1 / Constants.LAMBDA * (gradderivprime - gradderiv);
//
		System.err.println("finite " + gradderiv2test);
//
//
		System.err.println("----TESTING CONCLUDED----");


		return Utils.isSimilar(densityderiv2alpha, densityderiv2alphafinite, 1E-5) &&
				Utils.isSimilar(densityderiv2beta, densityderiv2betafinite, 1E-5) &&
				Math.abs(Hfderiv - Hfderivtest) < 1E-3 && Math.abs(IEderiv2 - IEtest) < 1E-3 &&
				Math.abs(dipolederiv2 - dipolederivtest) < 1E-3 && Math.abs(gradderiv2test - gradderiv2) < 1E-5;
	}

	public static boolean verifyEquations(SolutionU soln, int Z1, int param1) {


		SolutionU solnprime = (SolutionU) soln.withNewAtoms(Utils.perturbAtomParams(soln.atoms, Z1, param1));

		solnprime.compute();

		double Hfderiv = (solnprime.hf - soln.hf) / Constants.LAMBDA;

		double hfderivcheck = ParamDerivative.HfDeriv(soln, Z1, param1);
		System.err.println("---------");

		System.err.println(Hfderiv);

		System.err.println(hfderivcheck);
		System.err.println("---------");

		SimpleMatrix[][] matrices = ParamDerivative.staticDeriv(soln, Z1, 0);

		SimpleMatrix H = matrices[0][param1];

		SimpleMatrix Fa = matrices[1][param1];

		SimpleMatrix Fb = matrices[2][param1];

		System.err.println("mndohfderiv= " + ParamDerivative.HfDeriv(soln, H, Fa, Fb));

		SimpleMatrix xvector = PopleThiel.pt(soln, PopleThiel.toMO(soln.CtaOcc, soln.CaVirt, new SimpleMatrix[]{Fa}),
				PopleThiel.toMO(soln.CtbOcc, soln.CbVirt, new SimpleMatrix[]{Fb}))[0];

		SimpleMatrix[] densityderivs = PopleThiel.densityDeriv(soln, xvector);

		double dipolederiv = (solnprime.dipole - soln.dipole) / Constants.LAMBDA;

		double dipolederivcheck =
				ParamDerivative.nddoDipoleDeriv(soln, densityderivs[0].plus(densityderivs[1]), Z1, param1);
		System.err.println("---------");

		System.err.println(dipolederiv);

		System.err.println(dipolederivcheck);
		System.err.println("---------");

		SimpleMatrix[] responsematrices = PopleThiel.responseMatrix(soln, densityderivs);

		SimpleMatrix Fafull = soln.Cta.mult(Fa.plus(responsematrices[0])).mult(soln.Cta.transpose());

		SimpleMatrix Fbfull = soln.Ctb.mult(Fb.plus(responsematrices[1])).mult(soln.Ctb.transpose());

		SimpleMatrix[] x = ParamDerivative.xMatrix(soln, Fafull, Fbfull);

		SimpleMatrix Ca = soln.Cta.transpose().mult(x[0]);

		double homoderiv = (solnprime.homo - soln.homo) / Constants.LAMBDA;

		double homoderivcheck =
				ParamDerivative.homoDeriv(soln, x[0], soln.Cta.mult(Fa.plus(responsematrices[0])).mult(soln.Ca));
		System.err.println("---------");

		System.err.println(homoderiv);

		System.err.println(homoderivcheck);
		System.err.println("----expgeom start-----");


		double geomderiv =
				(GeometryDerivative.grad(solnprime, 0, 2) - GeometryDerivative.grad(soln, 0, 2)) / Constants.LAMBDA;

		double geomderivcheck =
				ParamGeometryDerivative.gradDeriv(soln, 0, 2, Z1, param1, densityderivs[0], densityderivs[1]);

		System.err.println(geomderiv);

		System.err.println(geomderivcheck);
		System.err.println("-----expgeom end----");

		return Math.abs(Hfderiv - hfderivcheck) < 1E-2 && Math.abs(homoderiv - homoderivcheck) < 1E-3 &&
				Math.abs(dipolederiv - dipolederivcheck) < 1E-3 && Math.abs(geomderiv - geomderivcheck) < 1E-5;
	}

	public static SimpleMatrix computeResponseVectorsThiel(SimpleMatrix xarray, SolutionU soln) {

		int NOccAlpha = soln.rm.nOccAlpha;
		int NOccBeta = soln.rm.nOccBeta;

		int NVirtAlpha = soln.rm.nVirtAlpha;
		int NVirtBeta = soln.rm.nVirtBeta;

		SimpleMatrix densityderivalpha = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		SimpleMatrix densityderivbeta = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);


		for (int u = 0; u < densityderivalpha.numRows(); u++) {
			for (int v = u; v < densityderivalpha.numCols(); v++) {
				double sum = 0;
				int count = 0;
				for (int i = 0; i < NOccAlpha; i++) {
					for (int j = 0; j < NVirtAlpha; j++) {
						sum -= (soln.Cta.get(i, u) * soln.Cta.get(j + NOccAlpha, v) +
								soln.Cta.get(j + NOccAlpha, u) * soln.Cta.get(i, v)) * xarray.get(count, 0);
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
								soln.Ctb.get(j + NOccBeta, u) * soln.Ctb.get(i, v)) * xarray.get(count, 0);
						count++;
					}
				}

				densityderivbeta.set(u, v, sum);
				densityderivbeta.set(v, u, sum);
			}
		}

		SimpleMatrix Jderiv = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		SimpleMatrix Kaderiv = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		SimpleMatrix Kbderiv = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);


		int Jcount = 0;
		int Kcount = 0;

		//construct J matrix

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double val = 0;
				if (j == k) {

					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						val += (densityderivalpha.get(l, l) + densityderivbeta.get(l, l)) *
								soln.integralArrayCoulomb[Jcount];
						Jcount++;
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += (densityderivalpha.get(l, m) + densityderivbeta.get(l, m)) *
										soln.integralArrayCoulomb[Jcount];
								Jcount++;
							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += (densityderivalpha.get(j, k) + densityderivbeta.get(j, k)) *
							soln.integralArrayCoulomb[Jcount];
					Jcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += (densityderivalpha.get(l, m) + densityderivbeta.get(l, m)) *
										soln.integralArrayCoulomb[Jcount];
								Jcount++;
							}
						}
					}
				}


				Jderiv.set(j, k, val);
				Jderiv.set(k, j, val);
			}
		}

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double vala = 0;
				double valb = 0;
				if (j == k) {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						vala += densityderivalpha.get(l, l) * soln.integralArrayExchange[Kcount];
						valb += densityderivbeta.get(l, l) * soln.integralArrayExchange[Kcount];
						Kcount++;
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					vala += densityderivalpha.get(j, k) * soln.integralArrayExchange[Kcount];
					valb += densityderivbeta.get(j, k) * soln.integralArrayExchange[Kcount];
					Kcount++;

				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.orbsOfAtom[soln.atomOfOrb[k]]) {
							vala += densityderivalpha.get(l, m) * soln.integralArrayExchange[Kcount];
							valb += densityderivbeta.get(l, m) * soln.integralArrayExchange[Kcount];
							Kcount++;
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


		SimpleMatrix R = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);

		int count1 = 0;

		for (int i = 0; i < NOccAlpha; i++) {
			for (int j = 0; j < NVirtAlpha; j++) {

				double element = 0;

				for (int u = 0; u < soln.nOrbitals; u++) {
					for (int v = 0; v < soln.nOrbitals; v++) {
						element += soln.Cta.get(i, u) * soln.Cta.get(j + NOccAlpha, v) * responsealpha.get(u, v);
					}
				}

				R.set(count1, 0, element);

				count1++;
			}
		}

		for (int i = 0; i < NOccBeta; i++) {
			for (int j = 0; j < NVirtBeta; j++) {

				double element = 0;

				for (int u = 0; u < soln.nOrbitals; u++) {
					for (int v = 0; v < soln.nOrbitals; v++) {
						element += soln.Ctb.get(i, u) * soln.Ctb.get(j + NOccBeta, v) * responsebeta.get(u, v);
					}
				}

				R.set(count1, 0, element);

				count1++;
			}
		}


		SimpleMatrix p = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);

		int counter = 0;

		for (int i = 0; i < NOccAlpha; i++) {
			for (int j = 0; j < NVirtAlpha; j++) {

				p.set(counter, 0,
						-R.get(counter, 0) + (soln.Ea.get(j + NOccAlpha) - soln.Ea.get(i)) * xarray.get(counter));
				counter++;
			}
		}

		for (int i = 0; i < NOccBeta; i++) {
			for (int j = 0; j < NVirtBeta; j++) {

				p.set(counter, 0,
						-R.get(counter, 0) + (soln.Eb.get(j + NOccBeta) - soln.Eb.get(i)) * xarray.get(counter));
				counter++;
			}
		}

		return p;
	}
}
