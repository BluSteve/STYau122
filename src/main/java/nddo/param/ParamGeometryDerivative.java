package nddo.param;

import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;

public class ParamGeometryDerivative {
	public static double gradderiv(SolutionR soln, int atomnum, int tau, int Z, int paramnum,
								   SimpleMatrix densityderiv) {
		double e = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			if (a != atomnum) {
				e += Ederiv(soln, atomnum, a, densityderiv, tau, Z, paramnum);
			}
		}

		return e;
	}

	public static double gradderiv(SolutionU soln, int atomnum, int tau, int Z, int paramnum,
								   SimpleMatrix densityderivalpha, SimpleMatrix densityderivbeta) {
		double e = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			if (a != atomnum) {
				e += Ederiv(soln, atomnum, a, densityderivalpha, densityderivbeta, tau, Z, paramnum);
			}
		}

		return e;
	}

	private static double Ederivstandard(SolutionR soln, int atomnum1, int atomnum2, SimpleMatrix densityderiv,
										 int tau) {

		int[][] index = soln.orbsOfAtom;

		SimpleMatrix densityMatrix = soln.densityMatrix();

		NDDOAtom[] atoms = soln.atoms;

		NDDOOrbital[] orbitals = soln.orbitals;

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += densityderiv.get(i, j) * atoms[atomnum2].Vgd(orbitals[i], orbitals[j], tau);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e -= densityderiv.get(k, l) * atoms[atomnum1].Vgd(orbitals[k], orbitals[l], tau);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * densityderiv.get(i, k) * State.nom.Hgd(orbitals[i], orbitals[k], tau);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {
							e += (densityderiv.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderiv.get(k, l) -
									densityderiv.get(i, k) * 0.5 * densityMatrix.get(j, l) -
									densityMatrix.get(i, k) * 0.5 * densityderiv.get(j, l))
									* State.nom.Ggd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
									tau);
						}
					}
				}
			}
		}

		return e;
	}

	private static double Ederivstandard(SolutionU soln, int atomnum1, int atomnum2, SimpleMatrix densityderivalpha,
										 SimpleMatrix densityderivbeta, int tau) {

		int[][] index = soln.orbsOfAtom;

		SimpleMatrix densityMatrix = soln.densityMatrix();

		SimpleMatrix densityderiv = densityderivalpha.plus(densityderivbeta);

		NDDOAtom[] atoms = soln.atoms;

		NDDOOrbital[] orbitals = soln.orbitals;

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += densityderiv.get(i, j) * atoms[atomnum2].Vgd(orbitals[i], orbitals[j], tau);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e -= densityderiv.get(k, l) * atoms[atomnum1].Vgd(orbitals[k], orbitals[l], tau);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * densityderiv.get(i, k) * State.nom.Hgd(orbitals[i], orbitals[k], tau);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {
							e += (densityderiv.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderiv.get(k, l) -
									densityderivalpha.get(i, k) * soln.alphaDensity().get(j, l) -
									densityderivbeta.get(i, k) * soln.betaDensity().get(j, l) -
									soln.alphaDensity().get(i, k) * densityderivalpha.get(j, l) -
									soln.betaDensity().get(i, k) * densityderivbeta.get(j, l))
									* State.nom.Ggd(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau);
						}
					}
				}
			}
		}

		return e;
	}

	private static double Ederiv(SolutionR soln, int atomnum1, int atomnum2, SimpleMatrix densityderiv, int tau, int Z,
								 int paramnum) {
		if (soln.atomicNumbers[atomnum1] != Z && soln.atomicNumbers[atomnum2] != Z) {
			return Ederivstandard(soln, atomnum1, atomnum2, densityderiv, tau);
		}

		else if (paramnum != 1 && paramnum != 2 && paramnum != 5 && paramnum != 6) {
			return Ederivstandard(soln, atomnum1, atomnum2, densityderiv, tau);
		}

		else if (paramnum == 1 || paramnum == 2) {//beta
			double e = Ederivstandard(soln, atomnum1, atomnum2, densityderiv, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					if (i != -1 && k != -1) {
						e += 2 * densityMatrix.get(i, k) *
								State.nom.Hbetapgd(orbitals[i], orbitals[k],
										ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
												soln.atomicNumbers[soln.atomOfOrb[k]], Z, orbitals[i].getL(),
												orbitals[k].getL(), paramnum - 1), tau);
					}
				}
			}

			return e;
		}

		else {

			double e = Ederivstandard(soln, atomnum1, atomnum2, densityderiv, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					if (i != -1 && j != -1) {
						e += densityMatrix.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
								ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z),
								paramnum - 5, tau);
					}
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					if (k != -1 && l != -1) {
						e -= densityMatrix.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
								ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z),
								paramnum - 5, tau);
					}
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					if (i != -1 && k != -1) {
						e += 2 * densityMatrix.get(i, k) * State.nom.Hzetapgd(orbitals[i], orbitals[k],
								ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z),
								paramnum - 5, tau);
					}
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							if (i != -1 && j != -1 && k != -1 && l != -1) {
								e += (densityMatrix.get(i, j) * densityMatrix.get(k, l) -
										densityMatrix.get(i, k) * 0.5 * densityMatrix.get(j, l))
										* State.nom.Gpgd(orbitals[i], orbitals[j], orbitals[k],
										orbitals[l], ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
												soln.atomicNumbers[atomnum2], Z), paramnum - 5, tau);
							}
						}
					}
				}
			}

			return e;
		}
	}

	private static double Ederiv(SolutionU soln, int atomnum1, int atomnum2, SimpleMatrix densityderivalpha,
								 SimpleMatrix densityderivbeta, int tau, int Z, int paramnum) {

		if (soln.atomicNumbers[atomnum1] != Z && soln.atomicNumbers[atomnum2] != Z) {
			return Ederivstandard(soln, atomnum1, atomnum2, densityderivalpha, densityderivbeta, tau);
		}

		else if (paramnum != 1 && paramnum != 2 && paramnum != 5 && paramnum != 6) {
			return Ederivstandard(soln, atomnum1, atomnum2, densityderivalpha, densityderivbeta, tau);
		}

		else if (paramnum == 1 || paramnum == 2) {//beta
			double e = Ederivstandard(soln, atomnum1, atomnum2, densityderivalpha, densityderivbeta, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					if (i != -1 && k != -1) {
						e += 2 * densityMatrix.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
								ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
										soln.atomicNumbers[soln.atomOfOrb[k]], Z, orbitals[i].getL(),
										orbitals[k].getL(), paramnum - 1), tau);
					}
				}
			}

			return e;
		}

		else {

			double e = Ederivstandard(soln, atomnum1, atomnum2, densityderivalpha, densityderivbeta, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					if (i != -1 && j != -1) {
						e += densityMatrix.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
								ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z),
								paramnum - 5, tau);
					}
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					if (k != -1 && l != -1) {
						e -= densityMatrix.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
								ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z),
								paramnum - 5, tau);
					}
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					if (i != -1 && k != -1) {
						e += 2 * densityMatrix.get(i, k) * State.nom.Hzetapgd(orbitals[i], orbitals[k],
								ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z),
								paramnum - 5, tau);
					}
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							if (i != -1 && j != -1 && k != -1 && l != -1) {
								e += (densityMatrix.get(i, j) * densityMatrix.get(k, l) -
										soln.alphaDensity().get(i, k) * soln.alphaDensity().get(j, l) -
										soln.betaDensity().get(i, k) * soln.betaDensity().get(j, l))
										* State.nom.Gpgd(orbitals[i], orbitals[j], orbitals[k],
										orbitals[l], ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
												soln.atomicNumbers[atomnum2], Z), paramnum - 5, tau);
							}
						}
					}
				}
			}

			return e;
		}
	}
}
