package nddo.param;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.geometry.GeometryDerivative;
import nddo.geometry.GeometrySecondDerivative;
import nddo.math.PopleThiel;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.apache.commons.lang3.time.StopWatch;
import org.ejml.data.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;
import tools.Pow;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;

public class ParamSecondDerivative {

	//NOTE: this class currently just verifies the accuracy of equations.

	public static SimpleMatrix xmatrix(SimpleMatrix F, SolutionR soln) {

		SimpleMatrix x = new SimpleMatrix(F.numRows(), F.numRows());

		for (int i = 0; i < F.numRows(); i++) {
			for (int j = i + 1; j < F.numRows(); j++) {
				x.set(i, j, F.get(i, j) / (soln.E.get(j) - soln.E.get(i)));
				x.set(j, i, F.get(i, j) / (soln.E.get(i) - soln.E.get(j)));
			}
		}

		return x;

	}

	public static SimpleMatrix[] xmatrices(SimpleMatrix Fa, SimpleMatrix Fb, SolutionU soln) {

		SimpleMatrix xa = new SimpleMatrix(Fa.numRows(), Fa.numRows());
		SimpleMatrix xb = new SimpleMatrix(Fa.numRows(), Fa.numRows());


		for (int i = 0; i < Fa.numRows(); i++) {
			for (int j = i + 1; j < Fa.numRows(); j++) {
				xa.set(i, j, Fa.get(i, j) / (soln.Ea.get(j) - soln.Ea.get(i)));
				xa.set(j, i, Fa.get(i, j) / (soln.Ea.get(i) - soln.Ea.get(j)));

				xb.set(i, j, Fb.get(i, j) / (soln.Eb.get(j) - soln.Eb.get(i)));
				xb.set(j, i, Fb.get(i, j) / (soln.Eb.get(i) - soln.Eb.get(j)));
			}
		}

		return new SimpleMatrix[]{xa, xb};

	}

	public static SimpleMatrix gammamatrix(SimpleMatrix totalderiv, SolutionR soln) {

		SimpleMatrix gamma = new SimpleMatrix(totalderiv.numRows(), totalderiv.numRows());

		for (int i = 0; i < totalderiv.numRows(); i++) {
			for (int j = i + 1; j < totalderiv.numRows(); j++) {
				gamma.set(i, j, totalderiv.get(i, j) / (soln.E.get(j) - soln.E.get(i)));
				gamma.set(j, i, totalderiv.get(j, i) / (soln.E.get(i) - soln.E.get(j)));
			}
		}

		return gamma;

	}

	public static SimpleMatrix gammamatrixalpha(SimpleMatrix totalderivalpha, SolutionU soln) {

		SimpleMatrix gamma = new SimpleMatrix(totalderivalpha.numRows(), totalderivalpha.numRows());

		for (int i = 0; i < totalderivalpha.numRows(); i++) {
			for (int j = i + 1; j < totalderivalpha.numRows(); j++) {
				gamma.set(i, j, totalderivalpha.get(i, j) / (soln.Ea.get(j) - soln.Ea.get(i)));
				gamma.set(j, i, totalderivalpha.get(j, i) / (soln.Ea.get(i) - soln.Ea.get(j)));
			}
		}

		return gamma;

	}


	public static SimpleMatrix densityderiv2static(SolutionR s, SimpleMatrix x1, SimpleMatrix x2) {

		int NOcc = (int) (s.nElectrons / 2.0);

		SimpleMatrix C = s.C;

		SimpleMatrix Dstatic = new SimpleMatrix(s.nOrbitals, s.nOrbitals);

		for (int u = 0; u < s.nOrbitals; u++) {
			for (int v = u; v < s.nOrbitals; v++) {
				double sum = 0;
				for (int i = 0; i < NOcc; i++) {
					for (int j = 0; j < s.nOrbitals; j++) {
						for (int k = 0; k < s.nOrbitals; k++) {
							sum += 2 * (C.get(u, k) * C.get(v, i) + C.get(u, i) * C.get(v, k)) *
									(x1.get(j, i) * x2.get(k, j) + x1.get(k, j) * x2.get(j, i));
							sum += 2 * (C.get(u, j) * C.get(v, k) + C.get(u, k) * C.get(v, j)) * x1.get(j, i) *
									x2.get(k, i);
						}
					}
				}

				for (int i = 0; i < NOcc; i++) {
					for (int j = 0; j < NOcc; j++) {
						for (int k = 0; k < s.nOrbitals; k++) {
							sum += 2 * C.get(u, i) * C.get(v, j) *
									(x1.get(k, i) * x2.get(k, j) + x1.get(k, j) * x2.get(k, i));
						}
					}
				}

				Dstatic.set(u, v, sum);

				Dstatic.set(v, u, sum);

			}
		}

		return Dstatic;


	}

	public static SimpleMatrix[] densityderiv2static(SolutionU soln, SimpleMatrix x1a, SimpleMatrix x1b,
													 SimpleMatrix x2a, SimpleMatrix x2b) {

		int NOccAlpha = soln.getRm().nOccAlpha;
		int NOccBeta = soln.getRm().nOccBeta;

		SimpleMatrix Ca = soln.ca.transpose();

		SimpleMatrix Dstaticalpha = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int u = 0; u < soln.nOrbitals; u++) {
			for (int v = u; v < soln.nOrbitals; v++) {
				double sum = 0;
				for (int i = 0; i < NOccAlpha; i++) {
					for (int j = 0; j < soln.nOrbitals; j++) {
						for (int k = 0; k < soln.nOrbitals; k++) {
							sum += (Ca.get(u, k) * Ca.get(v, i) + Ca.get(u, i) * Ca.get(v, k)) *
									(x1a.get(j, i) * x2a.get(k, j) + x1a.get(k, j) * x2a.get(j, i));
							sum += (Ca.get(u, j) * Ca.get(v, k) + Ca.get(u, k) * Ca.get(v, j)) * x1a.get(j, i) *
									x2a.get(k, i);
						}
					}
				}

				for (int i = 0; i < NOccAlpha; i++) {
					for (int j = 0; j < NOccAlpha; j++) {
						for (int k = 0; k < soln.nOrbitals; k++) {
							sum += Ca.get(u, i) * Ca.get(v, j) *
									(x1a.get(k, i) * x2a.get(k, j) + x1a.get(k, j) * x2a.get(k, i));
						}
					}
				}

				Dstaticalpha.set(u, v, sum);

				Dstaticalpha.set(v, u, sum);

			}
		}

		SimpleMatrix Cb = soln.cb.transpose();

		SimpleMatrix Dstaticbeta = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		for (int u = 0; u < soln.nOrbitals; u++) {
			for (int v = u; v < soln.nOrbitals; v++) {
				double sum = 0;
				for (int i = 0; i < NOccBeta; i++) {
					for (int j = 0; j < soln.nOrbitals; j++) {
						for (int k = 0; k < soln.nOrbitals; k++) {
							sum += (Cb.get(u, k) * Cb.get(v, i) + Cb.get(u, i) * Cb.get(v, k)) *
									(x1b.get(j, i) * x2b.get(k, j) + x1b.get(k, j) * x2b.get(j, i));
							sum += (Cb.get(u, j) * Cb.get(v, k) + Cb.get(u, k) * Cb.get(v, j)) * x1b.get(j, i) *
									x2b.get(k, i);
						}
					}
				}

				for (int i = 0; i < NOccBeta; i++) {
					for (int j = 0; j < NOccBeta; j++) {
						for (int k = 0; k < soln.nOrbitals; k++) {
							sum += Cb.get(u, i) * Cb.get(v, j) *
									(x1b.get(k, i) * x2b.get(k, j) + x1b.get(k, j) * x2b.get(k, i));
						}
					}
				}

				Dstaticbeta.set(u, v, sum);

				Dstaticbeta.set(v, u, sum);

			}
		}

		return new SimpleMatrix[]{Dstaticalpha, Dstaticbeta};


	}


	public static SimpleMatrix generalizedMNDOGderivstatic(SolutionR soln, SimpleMatrix mat, int Z, int paramNum) {

		if (paramNum < 5 || paramNum > 6) {
			return new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		}
		else {
			return zetaGderivstatic(soln, mat, Z, paramNum - 5);
		}
	}

	private static SimpleMatrix zetaGderivstatic(SolutionR soln, SimpleMatrix mat, int Z, int type) {

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
										sum += mat.get(l, m) *
												State.nom.Gpd(orbitals[j], orbitals[k], orbitals[l],
														orbitals[m],
														ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
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
									sum += mat.get(l, m) *
											(-0.5 * State.nom.Gpd(orbitals[j], orbitals[l], orbitals[k],
													orbitals[m],
													ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
															atomicnumbers[atomNumber[k]], Z), type));
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

	public static SimpleMatrix[] generalizedMNDOGderivstatic(SolutionU soln, SimpleMatrix alphamat,
															 SimpleMatrix betamat, int Z, int paramNum) {

		if (paramNum < 5 || paramNum > 6) {
			return new SimpleMatrix[]{new SimpleMatrix(soln.orbitals.length, soln.orbitals.length),
					new SimpleMatrix(soln.orbitals.length, soln.orbitals.length)};
		}
		else {
			return zetaGderivstatic(soln, alphamat, betamat, Z, paramNum - 5);
		}
	}

	private static SimpleMatrix[] zetaGderivstatic(SolutionU soln, SimpleMatrix alphamat, SimpleMatrix betamat, int Z,
												   int type) {


		NDDOOrbital[] orbitals = soln.orbitals;

		int[][] index = soln.orbsOfAtom;

		int[][] missingIndex = soln.missingOfAtom;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		SimpleMatrix Ga = new SimpleMatrix(orbitals.length, orbitals.length);

		SimpleMatrix Gb = new SimpleMatrix(orbitals.length, orbitals.length);

		SimpleMatrix mat = alphamat.plus(betamat);


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
										suma += mat.get(l, m) *
												State.nom.Gpd(orbitals[j], orbitals[k], orbitals[l],
														orbitals[m],
														ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z), type);
										sumb += mat.get(l, m) *
												State.nom.Gpd(orbitals[j], orbitals[k], orbitals[l],
														orbitals[m],
														ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
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
									suma -= alphamat.get(l, m) *
											State.nom.Gpd(orbitals[j], orbitals[l], orbitals[k],
													orbitals[m], ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
															atomicnumbers[atomNumber[k]], Z), type);

									sumb -= betamat.get(l, m) *
											State.nom.Gpd(orbitals[j], orbitals[l], orbitals[k],
													orbitals[m], ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
															atomicnumbers[atomNumber[k]], Z), type);

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

	//here are some observations regarding the H-matrix second derivatives.

	//observe that the following second derivative combinations have a zero H-matrix derivative:

	//(betas, betas), (betas, betap), (betas, Uss), (betas, Upp), (betap, Uss), (betap, Upp), (Uss, Upp), (Uss, zetas)
	// , (Uss, zetap), (Upp, zetas), (Upp, zetap)

	//the following second derivative combinations have nontrivial H-matrix derivatives:

	//(betas, zetas), (betas, zetap), (betap, zetas), (betap, zetap), (zetas, zetas), (zetas, zetap), (zetap, zetap)

	protected static SimpleMatrix betazetaHderiv2(Solution soln, int Z1, int type1, int Z2, int type2) {
		NDDOOrbital[] orbitals = soln.orbitals;

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] != atomNumber[k]) {
					double Huk = State.nom.Hbetazetap2d(orbitals[j], orbitals[k],
							ParamDerivative.getNumBeta(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]], Z1,
									orbitals[j].getL(), orbitals[k].getL(), type1),
							ParamDerivative.getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]], Z2),
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

		int[] atomNumber = soln.atomOfOrb;

		int[] atomicnumbers = soln.atomicNumbers;

		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (atomNumber[j] == atomNumber[k]) {
					double Huv = 0;

					for (int an = 0; an < atoms.length; an++) {
						if (atomNumber[j] != an) {
							Huv += atoms[an].Vp2d(orbitals[j], orbitals[k],
									ParamDerivative.getNum(atomicnumbers[an], atomicnumbers[atomNumber[j]], Z1), type1,
									ParamDerivative.getNum(atomicnumbers[an], atomicnumbers[atomNumber[j]], Z2),
									type2);
						}
					}
					H.set(j, k, Huv);
					H.set(k, j, Huv);
				}
				else { // case 3
					double Huk = State.nom.Hzetazetap2d(orbitals[j], orbitals[k],
							ParamDerivative.getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]], Z1),
							type1,
							ParamDerivative.getNum(atomicnumbers[atomNumber[j]], atomicnumbers[atomNumber[k]], Z2),
							type2);
					H.set(j, k, Huk);
					H.set(k, j, Huk);
				}
			}
		}

		return H;

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

	public static SimpleMatrix zetazetaGderiv2static(SolutionR soln, int Z1, int type1, int Z2, int type2) {


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
												State.nom.Gp2d(orbitals[j], orbitals[k], orbitals[l],
														orbitals[m],
														ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z1), type1,
														ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z2), type2);

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
									sum += soln.densityMatrix().get(l, m) * (-0.5 *
											State.nom.Gp2d(orbitals[j], orbitals[l], orbitals[k],
													orbitals[m], ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
															atomicnumbers[atomNumber[k]], Z1), type1,
													ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
															atomicnumbers[atomNumber[k]], Z2), type2));
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

	public static SimpleMatrix[] zetazetaGderiv2static(SolutionU soln, int Z1, int type1, int Z2, int type2) {


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
												State.nom.Gp2d(orbitals[j], orbitals[k], orbitals[l],
														orbitals[m],
														ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z1), type1,
														ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z2), type2);
										sumb += soln.densityMatrix().get(l, m) *
												State.nom.Gp2d(orbitals[j], orbitals[k], orbitals[l],
														orbitals[m],
														ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z1), type1,
														ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
																atomicnumbers[atomNumber[l]], Z2), type2);

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
									suma += soln.alphaDensity().get(l, m) *
											-State.nom.Gp2d(orbitals[j], orbitals[l], orbitals[k],
													orbitals[m], ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
															atomicnumbers[atomNumber[k]], Z1), type1,
													ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
															atomicnumbers[atomNumber[k]], Z2), type2);
									sumb += soln.betaDensity().get(l, m) *
											-State.nom.Gp2d(orbitals[j], orbitals[l], orbitals[k],
													orbitals[m], ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
															atomicnumbers[atomNumber[k]], Z1), type1,
													ParamDerivative.getNum(atomicnumbers[atomNumber[j]],
															atomicnumbers[atomNumber[k]], Z2), type2);

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

	public static SimpleMatrix Gderiv2static(SolutionR soln, int Z1, int param1, int Z2, int param2) {

		if (param1 >= 5 && param1 <= 6) {//zeta
			if (param2 >= 5 && param2 <= 6) {
				return ParamSecondDerivative.zetazetaGderiv2static(soln, Z1, param1 - 5, Z2, param2 - 5);
			}
		}
		return new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

	}

	public static SimpleMatrix[] Gderiv2static(SolutionU soln, int Z1, int param1, int Z2, int param2) {

		if (param1 >= 5 && param1 <= 6) {//zeta
			if (param2 >= 5 && param2 <= 6) {
				return ParamSecondDerivative.zetazetaGderiv2static(soln, Z1, param1 - 5, Z2, param2 - 5);
			}
		}
		return new SimpleMatrix[]{new SimpleMatrix(soln.nOrbitals, soln.nOrbitals),
				new SimpleMatrix(soln.nOrbitals, soln.nOrbitals)};

	}


	public static double MNDOHFDeriv(SolutionR soln, int Z1, int param1, int Z2, int param2, SimpleMatrix Hfirst,
									 SimpleMatrix Ffirst, SimpleMatrix densityderiv, int densityderivparamtype) {

		SimpleMatrix F = null;

		SimpleMatrix H = Hderiv2(soln, Z1, param1, Z2, param2);


		switch (densityderivparamtype) {
			case 0:
				F = H.plus(Gderiv2static(soln, Z1, param1, Z2, param2))
						.plus(generalizedMNDOGderivstatic(soln, densityderiv, Z1, param1));
				break;
			case 1:
				F = H.plus(Gderiv2static(soln, Z1, param1, Z2, param2))
						.plus(generalizedMNDOGderivstatic(soln, densityderiv, Z2, param2));
		}

		double e = 0;

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = 0; k < soln.orbitals.length; k++) {
				e += 0.5 * soln.densityMatrix().get(j, k) * (H.get(j, k) + F.get(j, k));
				e += 0.5 * densityderiv.get(j, k) * (Hfirst.get(j, k) + Ffirst.get(j, k));

			}
		}

		return e / 4.3363E-2;

	}

	public static double MNDOHFDeriv(SolutionU soln, int Z1, int param1, int Z2, int param2, SimpleMatrix Hfirst,
									 SimpleMatrix Fafirst, SimpleMatrix Fbfirst, SimpleMatrix densityderivalpha,
									 SimpleMatrix densityderivbeta, int densityderivparamtype) {

		SimpleMatrix Fa = null;

		SimpleMatrix Fb = null;

		SimpleMatrix H = Hderiv2(soln, Z1, param1, Z2, param2);


		switch (densityderivparamtype) {
			case 0:

				SimpleMatrix[] matrices1 = Gderiv2static(soln, Z1, param1, Z2, param2);

				SimpleMatrix[] matrices2 =
						generalizedMNDOGderivstatic(soln, densityderivalpha, densityderivbeta, Z1, param1);

				Fa = H.plus(matrices1[0]).plus(matrices2[0]);
				Fb = H.plus(matrices1[1]).plus(matrices2[1]);
				break;
			case 1:
				matrices1 = Gderiv2static(soln, Z1, param1, Z2, param2);

				matrices2 = generalizedMNDOGderivstatic(soln, densityderivalpha, densityderivbeta, Z2, param2);
				Fa = H.plus(matrices1[0]).plus(matrices2[0]);
				Fb = H.plus(matrices1[1]).plus(matrices2[1]);
		}

		double e = 0;

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = 0; k < soln.orbitals.length; k++) {
				e += 0.5 * soln.alphaDensity().get(j, k) * (H.get(j, k) + Fa.get(j, k));
				e += 0.5 * densityderivalpha.get(j, k) * (Hfirst.get(j, k) + Fafirst.get(j, k));

				e += 0.5 * soln.betaDensity().get(j, k) * (H.get(j, k) + Fb.get(j, k));
				e += 0.5 * densityderivbeta.get(j, k) * (Hfirst.get(j, k) + Fbfirst.get(j, k));

			}
		}

		return e / 4.3363E-2;

	}

	public static double MNDOIEDeriv2(SolutionR soln, SimpleMatrix xA, SimpleMatrix xB, SimpleMatrix totalderiv,
									  SimpleMatrix Fderiva, SimpleMatrix Fderivb, SimpleMatrix Fderiv2) {


		SimpleMatrix gammadiag = SimpleMatrix.diag(xA.mult(xB).plus(xB.mult(xA)).diag().getDDRM().data).scale(0.5);

		SimpleMatrix gammaremainder = gammamatrix(totalderiv, soln);

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

	public static double MNDOIEDeriv2(SolutionU soln, SimpleMatrix xAalpha, SimpleMatrix xBalpha,
									  SimpleMatrix totalderivalpha, SimpleMatrix Faderiva, SimpleMatrix Faderivb,
									  SimpleMatrix Faderiv2) {


		SimpleMatrix gammadiag =
				SimpleMatrix.diag(xAalpha.mult(xBalpha).plus(xBalpha.mult(xAalpha)).diag().getDDRM().data).scale(0.5);

		SimpleMatrix gammaremainder = gammamatrixalpha(totalderivalpha, soln);

		SimpleMatrix gamma = gammadiag.plus(gammaremainder);

		SimpleMatrix Fa = soln.Fa;

		SimpleMatrix Ca = soln.ca.transpose();

		SimpleMatrix Caderiva = Ca.mult(xAalpha);

		SimpleMatrix Caderivb = Ca.mult(xBalpha);

		SimpleMatrix Caderiv2 = Ca.mult(gamma);

		SimpleMatrix Ea = Caderiv2.transpose().mult(Fa.mult(Ca)).plus(Caderiva.transpose().mult(Faderivb.mult(Ca)))
				.plus(Caderiva.transpose().mult(Fa.mult(Caderivb)));

		Ea = Ea.plus(Caderivb.transpose().mult(Faderiva.mult(Ca))).plus(Ca.transpose().mult(Faderiv2.mult(Ca)))
				.plus(Ca.transpose().mult(Faderiva.mult(Caderivb)));

		Ea = Ea.plus(Caderivb.transpose().mult(Fa.mult(Caderiva))).plus(Ca.transpose().mult(Faderivb.mult(Caderiva)))
				.plus(Ca.transpose().mult(Fa.mult(Caderiv2)));

		return Ea.diag().get(soln.getRm().nOccAlpha - 1);
	}

	public static double MNDODipoleDeriv2(Solution soln, SimpleMatrix densityderiva, SimpleMatrix densityderivb,
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

		int[][] index = soln.orbsOfAtom;

		SimpleMatrix densityMatrix = soln.densityMatrix();

		double[] populationsderiva = new double[atoms.length];
		double[] populationsderivb = new double[atoms.length];
		double[] populationsderiv2 = new double[atoms.length];


		for (int j = 0; j < atoms.length; j++) {
			double suma = 0;
			double sumb = 0;
			double sum2 = 0;
			for (int k : index[j]) {
				if (k > -1) {
					suma += densityderiva.get(k, k);
					sumb += densityderivb.get(k, k);
					sum2 += densityderiv2.get(k, k);

				}
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
			chargedipa[0] += 2.5416 * populationsderiva[j] * (atoms[j].getCoordinates()[0] - com[0]);
			chargedipa[1] += 2.5416 * populationsderiva[j] * (atoms[j].getCoordinates()[1] - com[1]);
			chargedipa[2] += 2.5416 * populationsderiva[j] * (atoms[j].getCoordinates()[2] - com[2]);

			chargedipb[0] += 2.5416 * populationsderivb[j] * (atoms[j].getCoordinates()[0] - com[0]);
			chargedipb[1] += 2.5416 * populationsderivb[j] * (atoms[j].getCoordinates()[1] - com[1]);
			chargedipb[2] += 2.5416 * populationsderivb[j] * (atoms[j].getCoordinates()[2] - com[2]);

			chargedip2[0] += 2.5416 * populationsderiv2[j] * (atoms[j].getCoordinates()[0] - com[0]);
			chargedip2[1] += 2.5416 * populationsderiv2[j] * (atoms[j].getCoordinates()[1] - com[1]);
			chargedip2[2] += 2.5416 * populationsderiv2[j] * (atoms[j].getCoordinates()[2] - com[2]);
		}


		double[] hybriddipa = new double[]{0, 0, 0};
		double[] hybriddipb = new double[]{0, 0, 0};
		double[] hybriddip2 = new double[]{0, 0, 0};


		for (int j = 0; j < atoms.length; j++) {

			if (index[j].length > 1) {//exclude hydrogen
				hybriddipa[0] =
						hybriddipa[0] - 2.5416 * 2 * atoms[j].D1() * densityderiva.get(index[j][0], index[j][1]);
				hybriddipa[1] =
						hybriddipa[1] - 2.5416 * 2 * atoms[j].D1() * densityderiva.get(index[j][0], index[j][2]);
				hybriddipa[2] =
						hybriddipa[2] - 2.5416 * 2 * atoms[j].D1() * densityderiva.get(index[j][0], index[j][3]);

				if (atoms[j].getAtomProperties().getZ() == Z1) {
					hybriddipa[0] = hybriddipa[0] - 2.5416 * 2 * D1deriva * densityMatrix.get(index[j][0],
							index[j][1]);
					hybriddipa[1] = hybriddipa[1] - 2.5416 * 2 * D1deriva * densityMatrix.get(index[j][0],
							index[j][2]);
					hybriddipa[2] = hybriddipa[2] - 2.5416 * 2 * D1deriva * densityMatrix.get(index[j][0],
							index[j][3]);
				}

				hybriddipb[0] =
						hybriddipb[0] - 2.5416 * 2 * atoms[j].D1() * densityderivb.get(index[j][0], index[j][1]);
				hybriddipb[1] =
						hybriddipb[1] - 2.5416 * 2 * atoms[j].D1() * densityderivb.get(index[j][0], index[j][2]);
				hybriddipb[2] =
						hybriddipb[2] - 2.5416 * 2 * atoms[j].D1() * densityderivb.get(index[j][0], index[j][3]);

				if (atoms[j].getAtomProperties().getZ() == Z2) {
					hybriddipb[0] = hybriddipb[0] - 2.5416 * 2 * D1derivb * densityMatrix.get(index[j][0],
							index[j][1]);
					hybriddipb[1] = hybriddipb[1] - 2.5416 * 2 * D1derivb * densityMatrix.get(index[j][0],
							index[j][2]);
					hybriddipb[2] = hybriddipb[2] - 2.5416 * 2 * D1derivb * densityMatrix.get(index[j][0],
							index[j][3]);
				}

				hybriddip2[0] =
						hybriddip2[0] - 2.5416 * 2 * atoms[j].D1() * densityderiv2.get(index[j][0], index[j][1]);
				hybriddip2[1] =
						hybriddip2[1] - 2.5416 * 2 * atoms[j].D1() * densityderiv2.get(index[j][0], index[j][2]);
				hybriddip2[2] =
						hybriddip2[2] - 2.5416 * 2 * atoms[j].D1() * densityderiv2.get(index[j][0], index[j][3]);

				if (atoms[j].getAtomProperties().getZ() == Z1) {
					hybriddip2[0] = hybriddip2[0] - 2.5416 * 2 * D1deriva * densityderivb.get(index[j][0],
							index[j][1]);
					hybriddip2[1] = hybriddip2[1] - 2.5416 * 2 * D1deriva * densityderivb.get(index[j][0],
							index[j][2]);
					hybriddip2[2] = hybriddip2[2] - 2.5416 * 2 * D1deriva * densityderivb.get(index[j][0],
							index[j][3]);
				}

				if (atoms[j].getAtomProperties().getZ() == Z2) {
					hybriddip2[0] = hybriddip2[0] - 2.5416 * 2 * D1derivb * densityderiva.get(index[j][0],
							index[j][1]);
					hybriddip2[1] = hybriddip2[1] - 2.5416 * 2 * D1derivb * densityderiva.get(index[j][0],
							index[j][2]);
					hybriddip2[2] = hybriddip2[2] - 2.5416 * 2 * D1derivb * densityderiva.get(index[j][0],
							index[j][3]);
				}

				if (atoms[j].getAtomProperties().getZ() == Z1 && D1deriv2 != 0) {
					hybriddip2[0] = hybriddip2[0] - 2.5416 * 2 * D1deriv2 * densityMatrix.get(index[j][0],
							index[j][1]);
					hybriddip2[1] = hybriddip2[1] - 2.5416 * 2 * D1deriv2 * densityMatrix.get(index[j][0],
							index[j][2]);
					hybriddip2[2] = hybriddip2[2] - 2.5416 * 2 * D1deriv2 * densityMatrix.get(index[j][0],
							index[j][3]);
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
				soln.dipole
				- (dipoletota[0] * soln.dipoletot[0] + dipoletota[1] * soln.dipoletot[1] +
				dipoletota[2] * soln.dipoletot[2]) *
				(dipoletotb[0] * soln.dipoletot[0] + dipoletotb[1] * soln.dipoletot[1] +
						dipoletotb[2] * soln.dipoletot[2]) / Pow.pow(soln.dipole, 3);
	}

	public static SimpleMatrix staticMatrix(SolutionR soln, SimpleMatrix Fstatictotal, SimpleMatrix FstaticA,
											SimpleMatrix FstaticB, SimpleMatrix xvectorA, SimpleMatrix xvectorB,
											int Z1, int param1, int Z2, int param2) {

		SimpleMatrix densityA = ParamDerivative.densityDerivativeLimited(soln, xvectorA);
		SimpleMatrix densityB = ParamDerivative.densityDerivativeLimited(soln, xvectorB);

		SimpleMatrix GaB = ParamSecondDerivative.generalizedMNDOGderivstatic(soln, densityB, Z1, param1);
		SimpleMatrix GbA = ParamSecondDerivative.generalizedMNDOGderivstatic(soln, densityA, Z2, param2);

		SimpleMatrix FockA = FstaticA.plus(ParamDerivative.responseMatrix(soln, densityA));
		SimpleMatrix FockB = FstaticB.plus(ParamDerivative.responseMatrix(soln, densityB));

		SimpleMatrix FA = soln.Ct.mult(FockA.mult(soln.C));
		SimpleMatrix diagFA = SimpleMatrix.diag(FA.diag().getDDRM().data);
		SimpleMatrix FB = soln.Ct.mult(FockB.mult(soln.C));
		SimpleMatrix diagFB = SimpleMatrix.diag(FB.diag().getDDRM().data);


		SimpleMatrix xA = xmatrix(FA, soln);
		SimpleMatrix xB = xmatrix(FB, soln);

		SimpleMatrix omega =
				ParamDerivative.responseMatrix(soln, ParamSecondDerivative.densityderiv2static(soln, xA, xB));

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

	public static SimpleMatrix[] staticMatrix(SolutionU soln, SimpleMatrix Fstatictotalalpha,
											  SimpleMatrix Fstatictotalbeta, SimpleMatrix FstaticAalpha,
											  SimpleMatrix FstaticAbeta, SimpleMatrix FstaticBalpha,
											  SimpleMatrix FstaticBbeta, SimpleMatrix xvectorA, SimpleMatrix xvectorB,
											  int Z1, int param1, int Z2, int param2) {

		SimpleMatrix[] densitiesA = ParamDerivative.densityDerivatives(soln, xvectorA);

		SimpleMatrix densityAalpha = densitiesA[0];
		SimpleMatrix densityAbeta = densitiesA[1];

		SimpleMatrix[] densitiesB = ParamDerivative.densityDerivatives(soln, xvectorB);

		SimpleMatrix densityBalpha = densitiesB[0];
		SimpleMatrix densityBbeta = densitiesB[1];

		SimpleMatrix[] GaB =
				ParamSecondDerivative.generalizedMNDOGderivstatic(soln, densityBalpha, densityBbeta, Z1, param1);

		SimpleMatrix GaBalpha = GaB[0];

		SimpleMatrix GaBbeta = GaB[1];

		SimpleMatrix[] GbA =
				ParamSecondDerivative.generalizedMNDOGderivstatic(soln, densityAalpha, densityAbeta, Z2, param2);

		SimpleMatrix GbAalpha = GbA[0];

		SimpleMatrix GbAbeta = GbA[1];

		SimpleMatrix[] Ra = ParamDerivative.responseMatrices(soln, densitiesA);

		SimpleMatrix FockAalpha = FstaticAalpha.plus(Ra[0]);
		SimpleMatrix FockAbeta = FstaticAbeta.plus(Ra[1]);

		SimpleMatrix[] Rb = ParamDerivative.responseMatrices(soln, densitiesB);

		SimpleMatrix FockBalpha = FstaticBalpha.plus(Rb[0]);
		SimpleMatrix FockBbeta = FstaticBbeta.plus(Rb[1]);

		SimpleMatrix FAalpha = soln.ca.mult(FockAalpha).mult(soln.ca.transpose());
		SimpleMatrix FAbeta = soln.cb.mult(FockAbeta).mult(soln.cb.transpose());

		SimpleMatrix FBalpha = soln.ca.mult(FockBalpha).mult(soln.ca.transpose());
		SimpleMatrix FBbeta = soln.cb.mult(FockBbeta).mult(soln.cb.transpose());

		SimpleMatrix diagFAalpha = SimpleMatrix.diag(FAalpha.diag().getDDRM().data);
		SimpleMatrix diagFBalpha = SimpleMatrix.diag(FBalpha.diag().getDDRM().data);

		SimpleMatrix diagFAbeta = SimpleMatrix.diag(FAbeta.diag().getDDRM().data);
		SimpleMatrix diagFBbeta = SimpleMatrix.diag(FBbeta.diag().getDDRM().data);

		SimpleMatrix[] xA = xmatrices(FAalpha, FAbeta, soln);

		SimpleMatrix xAalpha = xA[0];
		SimpleMatrix xAbeta = xA[1];

		SimpleMatrix[] xB = xmatrices(FBalpha, FBbeta, soln);

		SimpleMatrix xBalpha = xB[0];
		SimpleMatrix xBbeta = xB[1];

		SimpleMatrix[] omega = ParamDerivative.responseMatrices(soln,
				ParamSecondDerivative.densityderiv2static(soln, xAalpha, xAbeta, xBalpha, xBbeta));

		SimpleMatrix omegaalpha = omega[0];

		SimpleMatrix omegabeta = omega[1];

		SimpleMatrix Phialpha = Fstatictotalalpha.plus(GaBalpha).plus(GbAalpha).plus(omegaalpha);
		SimpleMatrix Phibeta = Fstatictotalbeta.plus(GaBbeta).plus(GbAbeta).plus(omegabeta);

		SimpleMatrix matrixalpha = soln.ca.mult(Phialpha).mult(soln.ca.transpose());
		SimpleMatrix matrixbeta = soln.cb.mult(Phibeta).mult(soln.cb.transpose());

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


	public static SimpleMatrix staticFockDeriv(SolutionR soln, SimpleMatrix Fstatictotal, SimpleMatrix FstaticA,
											   SimpleMatrix FstaticB, SimpleMatrix xvectorA, SimpleMatrix xvectorB,
											   int Z1, int param1, int Z2, int param2) {

		SimpleMatrix densityA = ParamDerivative.densityDerivativeLimited(soln, xvectorA);
		SimpleMatrix densityB = ParamDerivative.densityDerivativeLimited(soln, xvectorB);

		SimpleMatrix GaB = ParamSecondDerivative.generalizedMNDOGderivstatic(soln, densityB, Z1, param1);
		SimpleMatrix GbA = ParamSecondDerivative.generalizedMNDOGderivstatic(soln, densityA, Z2, param2);

		SimpleMatrix FockA = FstaticA.plus(ParamDerivative.responseMatrix(soln, densityA));
		SimpleMatrix FockB = FstaticB.plus(ParamDerivative.responseMatrix(soln, densityB));

		SimpleMatrix FA = soln.Ct.mult(FockA.mult(soln.C));
		SimpleMatrix FB = soln.Ct.mult(FockB.mult(soln.C));


		SimpleMatrix xA = xmatrix(FA, soln);
		SimpleMatrix xB = xmatrix(FB, soln);

		SimpleMatrix omega =
				ParamDerivative.responseMatrix(soln, ParamSecondDerivative.densityderiv2static(soln, xA, xB));

		return Fstatictotal.plus(GaB).plus(GbA).plus(omega);

	}

	public static SimpleMatrix[] staticFockDeriv(SolutionU soln, SimpleMatrix Fstatictotalalpha,
												 SimpleMatrix Fstatictotalbeta, SimpleMatrix FstaticAalpha,
												 SimpleMatrix FstaticAbeta, SimpleMatrix FstaticBalpha,
												 SimpleMatrix FstaticBbeta, SimpleMatrix xvectorA,
												 SimpleMatrix xvectorB, int Z1, int param1, int Z2, int param2) {

		SimpleMatrix[] densitiesA = ParamDerivative.densityDerivatives(soln, xvectorA);

		SimpleMatrix densityAalpha = densitiesA[0];
		SimpleMatrix densityAbeta = densitiesA[1];

		SimpleMatrix[] densitiesB = ParamDerivative.densityDerivatives(soln, xvectorB);

		SimpleMatrix densityBalpha = densitiesB[0];
		SimpleMatrix densityBbeta = densitiesB[1];

		SimpleMatrix[] GaB =
				ParamSecondDerivative.generalizedMNDOGderivstatic(soln, densityBalpha, densityBbeta, Z1, param1);

		SimpleMatrix GaBalpha = GaB[0];

		SimpleMatrix GaBbeta = GaB[1];

		SimpleMatrix[] GbA =
				ParamSecondDerivative.generalizedMNDOGderivstatic(soln, densityAalpha, densityAbeta, Z2, param2);

		SimpleMatrix GbAalpha = GbA[0];

		SimpleMatrix GbAbeta = GbA[1];

		SimpleMatrix[] Ra = ParamDerivative.responseMatrices(soln, densitiesA);

		SimpleMatrix FockAalpha = FstaticAalpha.plus(Ra[0]);
		SimpleMatrix FockAbeta = FstaticAbeta.plus(Ra[1]);

		SimpleMatrix[] Rb = ParamDerivative.responseMatrices(soln, densitiesB);

		SimpleMatrix FockBalpha = FstaticBalpha.plus(Rb[0]);
		SimpleMatrix FockBbeta = FstaticBbeta.plus(Rb[1]);

		SimpleMatrix FAalpha = soln.ca.mult(FockAalpha).mult(soln.ca.transpose());
		SimpleMatrix FAbeta = soln.cb.mult(FockAbeta).mult(soln.cb.transpose());

		SimpleMatrix FBalpha = soln.ca.mult(FockBalpha).mult(soln.ca.transpose());
		SimpleMatrix FBbeta = soln.cb.mult(FockBbeta).mult(soln.cb.transpose());

		SimpleMatrix[] xA = xmatrices(FAalpha, FAbeta, soln);

		SimpleMatrix xAalpha = xA[0];
		SimpleMatrix xAbeta = xA[1];

		SimpleMatrix[] xB = xmatrices(FBalpha, FBbeta, soln);

		SimpleMatrix xBalpha = xB[0];
		SimpleMatrix xBbeta = xB[1];

		SimpleMatrix[] omega = ParamDerivative.responseMatrices(soln,
				ParamSecondDerivative.densityderiv2static(soln, xAalpha, xAbeta, xBalpha, xBbeta));

		SimpleMatrix omegaalpha = omega[0];
		SimpleMatrix omegabeta = omega[1];

		SimpleMatrix Phialpha = Fstatictotalalpha.plus(GaBalpha).plus(GbAalpha).plus(omegaalpha);
		SimpleMatrix Phibeta = Fstatictotalbeta.plus(GaBbeta).plus(GbAbeta).plus(omegabeta);

		return new SimpleMatrix[]{Phialpha, Phibeta};
	}

	public static SimpleMatrix utilitylazy(SolutionR soln, SimpleMatrix Fstatictotal, SimpleMatrix FstaticA,
										   SimpleMatrix FstaticB, SimpleMatrix xvectorA, SimpleMatrix xvectorB, int Z1,
										   int param1, int Z2, int param2) {

		SimpleMatrix densityA = ParamDerivative.densityDerivativeLimited(soln, xvectorA);
		SimpleMatrix densityB = ParamDerivative.densityDerivativeLimited(soln, xvectorB);

		SimpleMatrix FockA = FstaticA.plus(ParamDerivative.responseMatrix(soln, densityA));
		SimpleMatrix FockB = FstaticB.plus(ParamDerivative.responseMatrix(soln, densityB));

		SimpleMatrix FA = soln.Ct.mult(FockA.mult(soln.C));
		SimpleMatrix FB = soln.Ct.mult(FockB.mult(soln.C));

		SimpleMatrix xA = xmatrix(FA, soln);
		SimpleMatrix xB = xmatrix(FB, soln);

		SimpleMatrix densityderivstatic = ParamSecondDerivative.densityderiv2static(soln, xA, xB);

		return densityderivstatic;
	}

	public static SimpleMatrix[] utilitylazy(SolutionU soln, SimpleMatrix Fstatictotalalpha,
											 SimpleMatrix Fstatictotalbeta, SimpleMatrix FstaticAalpha,
											 SimpleMatrix FstaticAbeta, SimpleMatrix FstaticBalpha,
											 SimpleMatrix FstaticBbeta, SimpleMatrix xvectorA, SimpleMatrix xvectorB,
											 int Z1, int param1, int Z2, int param2) {

		SimpleMatrix[] densitiesA = ParamDerivative.densityDerivatives(soln, xvectorA);

		SimpleMatrix[] densitiesB = ParamDerivative.densityDerivatives(soln, xvectorB);

		SimpleMatrix[] Ra = ParamDerivative.responseMatrices(soln, densitiesA);

		SimpleMatrix FockAalpha = FstaticAalpha.plus(Ra[0]);
		SimpleMatrix FockAbeta = FstaticAbeta.plus(Ra[1]);

		SimpleMatrix[] Rb = ParamDerivative.responseMatrices(soln, densitiesB);

		SimpleMatrix FockBalpha = FstaticBalpha.plus(Rb[0]);
		SimpleMatrix FockBbeta = FstaticBbeta.plus(Rb[1]);

		SimpleMatrix FAalpha = soln.ca.mult(FockAalpha).mult(soln.ca.transpose());
		SimpleMatrix FAbeta = soln.cb.mult(FockAbeta).mult(soln.cb.transpose());

		SimpleMatrix FBalpha = soln.ca.mult(FockBalpha).mult(soln.ca.transpose());
		SimpleMatrix FBbeta = soln.cb.mult(FockBbeta).mult(soln.cb.transpose());


		SimpleMatrix[] xA = xmatrices(FAalpha, FAbeta, soln);

		SimpleMatrix xAalpha = xA[0];
		SimpleMatrix xAbeta = xA[1];

		SimpleMatrix[] xB = xmatrices(FBalpha, FBbeta, soln);

		SimpleMatrix xBalpha = xB[0];
		SimpleMatrix xBbeta = xB[1];


		return ParamSecondDerivative.densityderiv2static(soln, xAalpha, xAbeta, xBalpha, xBbeta);
	}

	public static SimpleMatrix[] gammaArrayThiel(SolutionU soln, SimpleMatrix[] matrixderivalpha,
												 SimpleMatrix[] matrixderivbeta) {

		StopWatch sw = new StopWatch();
		sw.start();


		int NOccAlpha = soln.getRm().nOccAlpha;
		int NOccBeta = soln.getRm().nOccBeta;

		int NVirtAlpha = soln.getRm().nVirtAlpha;
		int NVirtBeta = soln.getRm().nVirtBeta;

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
			SimpleMatrix[] densityderivs =
					new SimpleMatrix[matrixderivalpha.length];

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
					p.add(D.mult(GeometrySecondDerivative.computeResponseVectorsThiel(dirs[i], soln)));
				}


			}


			SimpleMatrix solver = new SimpleMatrix(p.size(), p.size());


			SimpleMatrix rhsvec = new SimpleMatrix(p.size(), rarray.length);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					SimpleMatrix rhs = new SimpleMatrix(p.size(), 1);

					for (int i = 0; i < rhs.numRows(); i++) {
						rhs.set(i, 0, 2 *
								rarray[a].transpose().mult(d.get(i)).get(0,
										0));

					}

					rhsvec.setColumn(a, 0, rhs.getDDRM().data);
				}
			}

			for (int i = 0; i < solver.numRows(); i++) {
				for (int j = i; j < solver.numRows(); j++) {

					double val = p.get(j).transpose().mult(d.get(i)).get(0,
							0) +
							p.get(i).transpose().mult(d.get(j)).get(0, 0);
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

	public static SimpleMatrix[] GammaArrayLimitedPople(SolutionR soln,
														SimpleMatrix[] staticMatrices) {

		int NOcc = (int) (soln.nElectrons / 2.0);
		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix[] gammaArray = new SimpleMatrix[staticMatrices.length];
		SimpleMatrix[] gammaArrayHold = new SimpleMatrix[staticMatrices.length];
		SimpleMatrix[] barray = new SimpleMatrix[staticMatrices.length];
		SimpleMatrix[] parray = new SimpleMatrix[staticMatrices.length];
		SimpleMatrix[] Farray = new SimpleMatrix[staticMatrices.length];
		SimpleMatrix[] rArray = new SimpleMatrix[staticMatrices.length];


		SimpleMatrix D = SimpleMatrix.identity(NOcc * NVirt);

		SimpleMatrix Dinv = SimpleMatrix.identity(NOcc * NVirt);

		for (int a = 0; a < gammaArray.length; a++) {
			SimpleMatrix F = new SimpleMatrix(NOcc * NVirt, 1);

			int count1 = 0;

			for (int i = 0; i < NOcc; i++) {
				for (int j = 0; j < NVirt; j++) {
					double element = staticMatrices[a].get(i, j + NOcc) / (soln.E.get(j + NOcc) - soln.E.get(i));

					F.set(count1, 0, element);
					count1++;
				}
			}

			F = D.mult(F);
			gammaArray[a] = new SimpleMatrix(NOcc * NVirt, 1);
			rArray[a] = gammaArray[a].copy();
			barray[a] = F.copy();
			Farray[a] = F.copy();
		}


		if (barray[0].numRows() == 0) {
			SimpleMatrix[] densityderivs =
					new SimpleMatrix[staticMatrices.length];

			for (int i = 0; i < densityderivs.length; i++) {
				densityderivs[i] = new SimpleMatrix(0, 0);
			}

			return densityderivs;
		}


		ArrayList<SimpleMatrix> b = new ArrayList<>();

		ArrayList<SimpleMatrix> p = new ArrayList<>();

		int[] iterable = new int[barray.length];


		SimpleMatrix F = new SimpleMatrix(NOcc * NVirt, Farray.length);
		for (int i = 0; i < Farray.length; i++) {
			F.setColumn(i, 0, Farray[i].getDDRM().data);
		}

		double[] oldrMags = new double[rArray.length];
		Arrays.fill(oldrMags, 10);
		SimpleMatrix responseMatrix = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		bigLoop:
		while (Utils.numIterable(iterable) > 0) {
			for (int i = 0; i < barray.length; i++) {
				for (int j = 0; j < i; j++) {
					barray[i] = barray[i].minus(barray[j].scale(barray[i].dot(barray[j]) / barray[j].dot(barray[j])));
				}
			}

			System.out.println("only " + Utils
					.numIterable(iterable) + " left to go!");

			for (int i = 0; i < barray.length; i++) {
				b.add(barray[i].copy());
				parray[i] = D.mult(
						PopleThiel.computeResponseVectorsPople(soln, Dinv.mult(barray[i].copy()), responseMatrix));
				p.add(parray[i].copy());
			}

			for (int i = 0; i < barray.length; i++) {
				SimpleMatrix newb = parray[i];

				for (int j = 0; j < b.size(); j++) {
					double num = b.get(j).transpose().mult(parray[i]).get(0) /
							b.get(j).transpose().mult(b.get(j)).get(0);
					newb = newb.minus(b.get(j).scale(num));
				}

				barray[i] = newb.copy();
			}

			SimpleMatrix B = new SimpleMatrix(NOcc * NVirt, b.size());
			SimpleMatrix P = new SimpleMatrix(NOcc * NVirt, b.size());

			for (int i = 0; i < b.size(); i++) {
				B.setColumn(i, 0, b.get(i).getDDRM().data);

				P.setColumn(i, 0, b.get(i).minus(p.get(i)).getDDRM().data);
			}


			SimpleMatrix lhs = B.transpose().mult(P);
			SimpleMatrix rhs = B.transpose().mult(F);
			SimpleMatrix alpha;

			try {
				alpha = lhs.solve(rhs);
			} catch (SingularMatrixException e) {
				alpha = SimpleMatrix.ones(lhs.numCols(), rhs.numCols());
			}

			for (int a = 0; a < gammaArray.length; a++) {
				rArray[a] = new SimpleMatrix(NOcc * NVirt, 1);
				gammaArray[a] = new SimpleMatrix(NOcc * NVirt, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {

					rArray[j] =
							rArray[j].plus(b.get(i).minus(p.get(i))
									.scale(alpha.get(i, j)));
					gammaArray[j] = gammaArray[j].plus(b.get(i).scale(alpha.get(i, j)));
				}
			}

			int xArrayHoldNN = Utils.numNotNull(gammaArrayHold);
			for (int j = 0; j < alpha.numCols(); j++) {
				rArray[j] = rArray[j].minus(Farray[j]);
				gammaArray[j] = Dinv.mult(gammaArray[j]);

				double mag = Utils.mag(rArray[j]);
				if (mag > oldrMags[j] || mag != mag) {
//					System.err.println("something has gone wrong");
//					System.out.println("xArrayHold = " +
//							Arrays.toString(xArrayHold));
//					System.out.println("xArray = " + Arrays.toString(xArray));
//					System.out.println("mag = " + mag);
//					System.out.println("oldrMags = " + oldrMags[j]);
					if (xArrayHoldNN == gammaArrayHold.length) {
						System.err.println(
								"Some numerical instability encountered; " +
										"returning lower precision values...");
						gammaArray = gammaArrayHold;
						break bigLoop;
					}
					else {
						if (mag > oldrMags[j]) {
							System.err.println(
									"Numerical instability detected; " +
											"reverting to Thiel algorithm (which doesn't exist)...");
							System.exit(0);
						}
						if (mag != mag) {
							System.err.println("Pople algorithm fails; " +
									"reverting to Thiel algorithm (which doesn't exist)...");
							System.exit(0);
						}
					}
				}
				else {
					if (mag < 1E-7) {
						gammaArrayHold[j] = gammaArray[j];
						if (mag < 1E-8) {
							iterable[j] = 1;
						}
					}
					else {
						iterable[j] = 0;

					}
					System.out.println("Pople convergence test: " + mag);
				}

				oldrMags[j] = mag;
			}
		}
		return gammaArray;
	}


	public static boolean verifyEquations(SolutionR soln, int Z1, int param1, int Z2, int param2) {

		SolutionR solnprime = (SolutionR) soln.withNewAtoms(Utils.perturbAtomParams(soln.atoms, Z1, param1));

		solnprime.compute();

		SimpleMatrix H1 = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z1, 0)[0][param1];
		SimpleMatrix H2 = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z2, 0)[0][param2];


		SimpleMatrix F2 = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z2, 0)[1][param2];
		SimpleMatrix F2prime = ParamDerivative.MNDOStaticMatrixDeriv(solnprime, Z2, 0)[1][param2];
		SimpleMatrix F1 = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z1, 0)[1][param1];

		SimpleMatrix x2 = PopleThiel.pople(soln, new SimpleMatrix[]{F2})[0];
		SimpleMatrix x2prime = PopleThiel.pople(solnprime, new SimpleMatrix[]{F2prime})[0];
		SimpleMatrix x1 = PopleThiel.pople(soln, new SimpleMatrix[]{F1})[0];

		SimpleMatrix D1 = ParamDerivative.densityDerivativeLimited(soln, x1);

		SimpleMatrix R1 = ParamDerivative.responseMatrix(soln, D1);

		SimpleMatrix C1 = soln.C.mult(xmatrix(soln.Ct.mult(F1.plus(R1)).mult(soln.C), soln));


		SimpleMatrix D2prime = ParamDerivative.densityDerivativeLimited(solnprime, x2prime);


		SimpleMatrix D2 = ParamDerivative.densityDerivativeLimited(soln, x2);

		SimpleMatrix R2 = ParamDerivative.responseMatrix(soln, D2);

		SimpleMatrix C2 = soln.C.mult(xmatrix(soln.Ct.mult(F2.plus(R2)).mult(soln.C), soln));

		SimpleMatrix densityderiv2finite = D2prime.minus(D2).scale(1 / Constants.LAMBDA);

		SimpleMatrix Fstatictotal = ParamSecondDerivative.Hderiv2(soln, Z1, param1, Z2, param2)
				.plus(ParamSecondDerivative.Gderiv2static(soln, Z1, param1, Z2, param2));

		SimpleMatrix Hstatic = ParamSecondDerivative.Hderiv2(soln, Z1, param1, Z2, param2);


		SimpleMatrix densityderiv2static =
				ParamSecondDerivative.utilitylazy(soln, Fstatictotal, F1, F2, x1, x2, Z1, param1, Z2, param2);

		SimpleMatrix rhsmat =
				ParamSecondDerivative.staticMatrix(soln, Fstatictotal, F1, F2, x1, x2, Z1, param1, Z2, param2);

		SimpleMatrix gammavec = GammaArrayLimitedPople(soln, new SimpleMatrix[]{rhsmat})[0];

		SimpleMatrix densityderiv2response = ParamDerivative.densityDerivativeLimited(soln, gammavec);

		SimpleMatrix densityderiv2 = densityderiv2response.plus(densityderiv2static);

		double Hfderiv = 1 / Constants.LAMBDA *
				(ParamDerivative.HFDeriv(solnprime, Z2, param2) - ParamDerivative.HFDeriv(soln, Z2, param2));

		double Hfderivtest = ParamSecondDerivative.MNDOHFDeriv(soln, Z1, param1, Z2, param2, H1, F1, D2, 0);

		double Hfderivtest2 = ParamSecondDerivative.MNDOHFDeriv(soln, Z1, param1, Z2, param2, H2, F2, D1, 1);

		SimpleMatrix x1mat =
				xmatrix(soln.Ct.mult(F1.plus(ParamDerivative.responseMatrix(soln, D1))).mult(soln.C), soln);

		SimpleMatrix x2mat =
				xmatrix(soln.Ct.mult(F2.plus(ParamDerivative.responseMatrix(soln, D2))).mult(soln.C), soln);

		SimpleMatrix x2matprime =
				xmatrix(solnprime.Ct.mult(F2prime.plus(ParamDerivative.responseMatrix(solnprime, D2prime))).mult(solnprime.C), solnprime);

		SimpleMatrix totalderiv = rhsmat.plus(
				soln.Ct.mult(ParamDerivative.responseMatrix(soln, densityderiv2response)).mult(soln.C));

		SimpleMatrix Fderiv2 = staticFockDeriv(soln, Fstatictotal, F1, F2, x1, x2, Z1, param1, Z2, param2).plus(
				ParamDerivative.responseMatrix(soln, densityderiv2response));

		SimpleMatrix Fderiva = F1.plus(ParamDerivative.responseMatrix(soln, D1));

		SimpleMatrix Fderivb = F2.plus(ParamDerivative.responseMatrix(soln, D2));

		SimpleMatrix Fderivbprime = F2prime.plus(ParamDerivative.responseMatrix(solnprime, D2prime));

		double IEtest = MNDOIEDeriv2(soln, x1mat, x2mat, totalderiv, Fderiva, Fderivb, Fderiv2);

		double IEderiv = ParamDerivative.MNDOHomoDerivNew(soln, x2mat, Fderivb);

		double IEderivprime = ParamDerivative.MNDOHomoDerivNew(solnprime, x2matprime, Fderivbprime);


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

		double dipolederivtest = MNDODipoleDeriv2(soln, D1, D2, densityderiv2, Z1, param1, Z2, param2);

		double dipolederiv = ParamDerivative.MNDODipoleDeriv(soln, D2, Z2, param2);

		double dipolederivprime = ParamDerivative.MNDODipoleDeriv(solnprime, D2prime, Z2, param2);

		double dipolederiv2 = 1 / Constants.LAMBDA * (dipolederivprime - dipolederiv);

		System.err.println("----TESTING DIPOLE DERIVATIVES----");


		System.err.println(dipolederiv2);
		System.err.println(dipolederivtest);


		System.err.println("----TESTING EXPGEOM DERIVATIVES----");

		double gradderiv = ParamGeometryDerivative.gradderiv(soln, 0, 0, Z2, param2, D2);

		double gradderivprime = ParamGeometryDerivative.gradderiv(solnprime, 0, 0, Z2, param2, D2prime);

		double gradderiv2 =
				ParamGeometrySecondDerivative.gradderiv2(soln, 0, 0, Z1, param1, Z2, param2, D1, D2, densityderiv2);

		System.err.println(gradderiv2);

		double gradderiv2test = 1 / Constants.LAMBDA * (gradderivprime - gradderiv);

		System.err.println(gradderiv2test);


		System.err.println("----TESTING CONCLUDED----");


		return Solution.isSimilar(densityderiv2, densityderiv2finite, 1E-5) && Math.abs(Hfderiv - Hfderivtest) < 1E-3 &&
				Math.abs(IEderiv2 - IEtest) < 1E-3 && Math.abs(dipolederiv2 - dipolederivtest) < 1E-3 &&
				Math.abs(gradderiv2test - gradderiv2) < 1E-5;
	}

	public static boolean verifyEquations(SolutionU soln, int Z1, int param1, int Z2, int param2) {

		SolutionU solnprime = (SolutionU) soln.withNewAtoms(Utils.perturbAtomParams(soln.atoms, Z1, param1));

		solnprime.compute();

		SimpleMatrix H1 = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z1, 0)[0][param1];
		SimpleMatrix H2 = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z2, 0)[0][param2];


		SimpleMatrix F2alpha = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z2, 0)[1][param2];
		SimpleMatrix F2beta = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z2, 0)[2][param2];

		SimpleMatrix F2alphaprime = ParamDerivative.MNDOStaticMatrixDeriv(solnprime, Z2, 0)[1][param2];
		SimpleMatrix F2betaprime = ParamDerivative.MNDOStaticMatrixDeriv(solnprime, Z2, 0)[2][param2];

		SimpleMatrix F1alpha = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z1, 0)[1][param1];
		SimpleMatrix F1beta = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z1, 0)[2][param1];


		SimpleMatrix x2 =
				ParamDerivative.xArrayThiel(soln, new SimpleMatrix[]{F2alpha}, new SimpleMatrix[]{F2beta})[0];
		SimpleMatrix x2prime = ParamDerivative.xArrayThiel(solnprime, new SimpleMatrix[]{F2alphaprime},
				new SimpleMatrix[]{F2betaprime})[0];
		SimpleMatrix x1 =
				ParamDerivative.xArrayThiel(soln, new SimpleMatrix[]{F1alpha}, new SimpleMatrix[]{F1beta})[0];

		SimpleMatrix[] D1 = ParamDerivative.densityDerivatives(soln, x1);

		SimpleMatrix[] R1 = ParamDerivative.responseMatrices(soln, D1);

		SimpleMatrix[] x1matrix = xmatrices(soln.ca.mult(F1alpha.plus(R1[0])).mult(soln.ca.transpose()),
				soln.cb.mult(F1beta.plus(R1[1])).mult(soln.cb.transpose()), soln);

		SimpleMatrix C1alpha = soln.ca.transpose().mult(x1matrix[0]);
		SimpleMatrix C1beta = soln.cb.transpose().mult(x1matrix[1]);

		SimpleMatrix[] D2prime = ParamDerivative.densityDerivatives(solnprime, x2prime);

		SimpleMatrix[] D2 = ParamDerivative.densityDerivatives(soln, x2);

		SimpleMatrix[] R2 = ParamDerivative.responseMatrices(soln, D2);

		SimpleMatrix[] R2prime = ParamDerivative.responseMatrices(solnprime, D2prime);


		SimpleMatrix[] x2matrix = xmatrices(soln.ca.mult(F2alpha.plus(R2[0])).mult(soln.ca.transpose()),
				soln.cb.mult(F2beta.plus(R2[1])).mult(soln.cb.transpose()), soln);

		SimpleMatrix[] x2primematrix =
				xmatrices(solnprime.ca.mult(F2alphaprime.plus(R2prime[0])).mult(solnprime.ca.transpose()),
						solnprime.cb.mult(F2betaprime.plus(R2prime[1])).mult(solnprime.cb.transpose()), solnprime);


		SimpleMatrix C2alpha = soln.ca.transpose().mult(x2matrix[0]);
		SimpleMatrix C2beta = soln.cb.transpose().mult(x2matrix[1]);

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

		SimpleMatrix gammavec = gammaArrayThiel(soln, new SimpleMatrix[]{rhsmat[0]}, new SimpleMatrix[]{rhsmat[1]})[0];

		SimpleMatrix[] densityderiv2response = ParamDerivative.densityDerivatives(soln, gammavec);

		SimpleMatrix densityderiv2alpha = densityderiv2response[0].plus(densityderiv2static[0]);

		SimpleMatrix densityderiv2beta = densityderiv2response[1].plus(densityderiv2static[1]);

		System.err.println(densityderiv2alpha);

		System.err.println(densityderiv2alphafinite);

		System.err.println("---");

		System.err.println(densityderiv2beta);

		System.err.println(densityderiv2betafinite);

		double Hfderiv = 1 / Constants.LAMBDA *
				(ParamDerivative.HFDeriv(solnprime, Z2, param2) - ParamDerivative.HFDeriv(soln, Z2, param2));

		double Hfderivtest =
				ParamSecondDerivative.MNDOHFDeriv(soln, Z1, param1, Z2, param2, H1, F1alpha, F1beta, D2[0], D2[1], 0);

		double Hfderivtest2 =
				ParamSecondDerivative.MNDOHFDeriv(soln, Z1, param1, Z2, param2, H2, F2alpha, F2beta, D1[0], D1[1], 1);

		SimpleMatrix totalderiv = rhsmat[0].plus(
				soln.ca.mult(ParamDerivative.responseMatrices(soln, densityderiv2response)[0])
						.mult(soln.ca.transpose()));

		SimpleMatrix Fderiv2 =
				staticFockDeriv(soln, Fstatictotalalpha, Fstatictotalbeta, F1alpha, F1beta, F2alpha, F2beta, x1, x2,
						Z1,
						param1, Z2, param2)[0].plus(ParamDerivative.responseMatrices(soln, densityderiv2response)[0]);

		SimpleMatrix Fderiva = F1alpha.plus(R1[0]);

		SimpleMatrix Fderivb = F2alpha.plus(R2[0]);

		SimpleMatrix Fderivbprime = F2alphaprime.plus(ParamDerivative.responseMatrices(solnprime, D2prime)[0]);
//
		double IEtest = MNDOIEDeriv2(soln, x1matrix[0], x2matrix[0], totalderiv, Fderiva, Fderivb, Fderiv2);



		double IEderiv = ParamDerivative.MNDOHomoDerivNew(soln, x2matrix[0], Fderivb);
		double IEderivprime =
				ParamDerivative.MNDOHomoDerivNew(solnprime, x2primematrix[0], Fderivbprime);
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
				MNDODipoleDeriv2(soln, D1[0].plus(D1[1]), D2[0].plus(D2[1]),
						densityderiv2alpha.plus(densityderiv2beta),
						Z1, param1, Z2, param2);
//
		double dipolederiv = ParamDerivative.MNDODipoleDeriv(soln, D2[0].plus(D2[1]), Z2, param2);
//
		double dipolederivprime = ParamDerivative.MNDODipoleDeriv(solnprime, D2prime[0].plus(D2prime[1]), Z2, param2);
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
		double gradderiv = ParamGeometryDerivative.gradderiv(soln, 0, 0, Z2, param2, D2[0], D2[1]);
//
		double gradderivprime =
				ParamGeometryDerivative.gradderiv(solnprime, 0, 0, Z2, param2, D2prime[0], D2prime[1]);
//
		double gradderiv2 =
				ParamGeometrySecondDerivative.gradderiv2(soln, 0, 0, Z1, param1, Z2, param2, D1[0], D1[1], D2[0],
						D2[1], densityderiv2alpha, densityderiv2beta);
//
		System.err.println("analyt " + gradderiv2);
//
		double gradderiv2test = 1 / Constants.LAMBDA * (gradderivprime - gradderiv);
//
		System.err.println("finite " + gradderiv2test);
//
//
		System.err.println("----TESTING CONCLUDED----");


		return Solution.isSimilar(densityderiv2alpha, densityderiv2alphafinite, 1E-5) &&
				Solution.isSimilar(densityderiv2beta, densityderiv2betafinite, 1E-5) &&
				Math.abs(Hfderiv - Hfderivtest) < 1E-3 && Math.abs(IEderiv2 - IEtest) < 1E-3 &&
				Math.abs(dipolederiv2 - dipolederivtest) < 1E-3 &&Math.abs(gradderiv2test - gradderiv2) < 1E-5;
	}

	public static boolean verifyEquations(SolutionU soln, int Z1, int param1) {

		SolutionU solnprime = (SolutionU) soln.withNewAtoms(Utils.perturbAtomParams(soln.atoms, Z1, param1));

		solnprime.compute();

		double Hfderiv = (solnprime.hf - soln.hf) / Constants.LAMBDA;

		double hfderivcheck = ParamDerivative.HFDeriv(soln, Z1, param1);
		System.err.println ("---------");

		System.err.println(Hfderiv);

		System.err.println(hfderivcheck);
		System.err.println ("---------");

		SimpleMatrix[][] matrices = ParamDerivative.MNDOStaticMatrixDeriv(soln, Z1, 0);

		SimpleMatrix H = matrices[0][param1];

		SimpleMatrix Fa = matrices[1][param1];

		SimpleMatrix Fb = matrices[2][param1];

		SimpleMatrix xvector = ParamDerivative.xArrayThiel(soln, new SimpleMatrix[]{Fa}, new SimpleMatrix[]{Fb})[0];

		SimpleMatrix[] densityderivs = ParamDerivative.densityDerivatives(soln, xvector);

		double dipolederiv = (solnprime.dipole - soln.dipole) / Constants.LAMBDA;

		double dipolederivcheck =
				ParamDerivative.MNDODipoleDeriv(soln, densityderivs[0].plus(densityderivs[1]), Z1, param1);
		System.err.println ("---------");

		System.err.println(dipolederiv);

		System.err.println(dipolederivcheck);
		System.err.println ("---------");

		SimpleMatrix[] responsematrices = ParamDerivative.responseMatrices(soln, densityderivs);

		SimpleMatrix Fafull = soln.ca.mult(Fa.plus(responsematrices[0])).mult(soln.ca.transpose());

		SimpleMatrix Fbfull = soln.cb.mult(Fb.plus(responsematrices[1])).mult(soln.cb.transpose());

		SimpleMatrix[] x = xmatrices(Fafull, Fbfull, soln);

		SimpleMatrix Ca = soln.ca.transpose().mult(x[0]);

		double homoderiv = (solnprime.homo - soln.homo) / Constants.LAMBDA;

		double homoderivcheck = ParamDerivative.MNDOHomoDerivNew(soln, x[0], Fa.plus(responsematrices[0]));
		System.err.println ("---------");

		System.err.println(homoderiv);

		System.err.println(homoderivcheck);
		System.err.println ("----expgeom start-----");


		double geomderiv = (GeometryDerivative.grad(solnprime, 0, 2) - GeometryDerivative.grad(soln, 0, 2)) /
		Constants.LAMBDA;

		double geomderivcheck = ParamGeometryDerivative.gradderiv(soln, 0, 2, Z1, param1, densityderivs[0],
				densityderivs[1]);

		System.err.println (geomderiv);

		System.err.println (geomderivcheck);
		System.err.println ("-----expgeom end----");


		return Math.abs(Hfderiv - hfderivcheck) < 1E-2;
	}
}
