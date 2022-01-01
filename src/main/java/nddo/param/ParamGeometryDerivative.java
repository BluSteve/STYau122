package nddo.param;

import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;

public class ParamGeometryDerivative {
	public static double gradDerivAlpha(Solution soln, int atomNum, int tau, int Z) {
		double e = 0;
		
		for (int a = 0; a < soln.atoms.length; a++) {
			if (a != atomNum) {
				e += soln.atoms[atomNum].crfalphapgd(soln.atoms[a],
						ParamDerivative.getNum(soln.atoms[atomNum].getAtomProperties().getZ(),
								soln.atoms[a].getAtomProperties().getZ(), Z),
						tau);
			}
		}
		
		return e;
	}

	public static double gradDeriv(SolutionR soln, int atomNum, int tau, int Z, int paramNum,
								   SimpleMatrix densityDeriv) {
		double e = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			if (a != atomNum) {
				e += Ederiv(soln, atomNum, a, densityDeriv, tau, Z, paramNum);
			}
		}

		return e;
	}

	public static double gradDeriv(SolutionU soln, int atomNum, int tau, int Z, int paramNum,
								   SimpleMatrix densityDerivAlpha, SimpleMatrix densityDerivBeta) {
		double e = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			if (a != atomNum) {
				e += Ederiv(soln, atomNum, a, densityDerivAlpha, densityDerivBeta, tau, Z, paramNum);
			}
		}

		return e;
	}

	private static double Ederivstandard(SolutionR soln, int atomNum1, int atomNum2, SimpleMatrix densityDeriv,
										 int tau) {
		int[][] orbsOfAtom = soln.orbsOfAtom;
		SimpleMatrix densityMatrix = soln.densityMatrix();
		NDDOAtom[] atoms = soln.atoms;
		NDDOOrbital[] orbitals = soln.orbitals;

		double e = 0;

		for (int i : orbsOfAtom[atomNum1]) {
			for (int j : orbsOfAtom[atomNum1]) {
				e += densityDeriv.get(i, j) * atoms[atomNum2].Vgd(orbitals[i], orbitals[j], tau);
			}
		}

		for (int k : orbsOfAtom[atomNum2]) {
			for (int l : orbsOfAtom[atomNum2]) {
				e -= densityDeriv.get(k, l) * atoms[atomNum1].Vgd(orbitals[k], orbitals[l], tau);
			}
		}

		for (int i : orbsOfAtom[atomNum1]) {
			for (int k : orbsOfAtom[atomNum2]) {
				e += 2 * densityDeriv.get(i, k) * State.nom.Hgd(orbitals[i], orbitals[k], tau);
			}
		}

		for (int i : orbsOfAtom[atomNum1]) {
			for (int j : orbsOfAtom[atomNum1]) {
				for (int k : orbsOfAtom[atomNum2]) {
					for (int l : orbsOfAtom[atomNum2]) {
						e += (densityDeriv.get(i, j) * densityMatrix.get(k, l) +
								densityMatrix.get(i, j) * densityDeriv.get(k, l) -
								densityDeriv.get(i, k) * 0.5 * densityMatrix.get(j, l) -
								densityMatrix.get(i, k) * 0.5 * densityDeriv.get(j, l))
								* State.nom.Ggd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
								tau);
					}
				}
			}
		}

		return e;
	}

	private static double Ederivstandard(SolutionU soln, int atomNum1, int atomNum2, SimpleMatrix densityDerivAlpha,
										 SimpleMatrix densityDerivBeta, int tau) {

		int[][] orbsOfAtom = soln.orbsOfAtom;
		SimpleMatrix densityMatrix = soln.densityMatrix();
		SimpleMatrix densityDeriv = densityDerivAlpha.plus(densityDerivBeta);
		NDDOAtom[] atoms = soln.atoms;
		NDDOOrbital[] orbitals = soln.orbitals;

		double e = 0;

		for (int i : orbsOfAtom[atomNum1]) {
			for (int j : orbsOfAtom[atomNum1]) {
				e += densityDeriv.get(i, j) * atoms[atomNum2].Vgd(orbitals[i], orbitals[j], tau);
			}
		}

		for (int k : orbsOfAtom[atomNum2]) {
			for (int l : orbsOfAtom[atomNum2]) {
				e -= densityDeriv.get(k, l) * atoms[atomNum1].Vgd(orbitals[k], orbitals[l], tau);
			}
		}

		for (int i : orbsOfAtom[atomNum1]) {
			for (int k : orbsOfAtom[atomNum2]) {
				e += 2 * densityDeriv.get(i, k) * State.nom.Hgd(orbitals[i], orbitals[k], tau);
			}
		}

		for (int i : orbsOfAtom[atomNum1]) {
			for (int j : orbsOfAtom[atomNum1]) {
				for (int k : orbsOfAtom[atomNum2]) {
					for (int l : orbsOfAtom[atomNum2]) {
						e += (densityDeriv.get(i, j) * densityMatrix.get(k, l) +
								densityMatrix.get(i, j) * densityDeriv.get(k, l) -
								densityDerivAlpha.get(i, k) * soln.alphaDensity().get(j, l) -
								densityDerivBeta.get(i, k) * soln.betaDensity().get(j, l) -
								soln.alphaDensity().get(i, k) * densityDerivAlpha.get(j, l) -
								soln.betaDensity().get(i, k) * densityDerivBeta.get(j, l))
								* State.nom.Ggd(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau);
					}
				}
			}
		}

		return e;
	}

	private static double Ederiv(SolutionR soln, int atomNum1, int atomNum2, SimpleMatrix densityDeriv, int tau, int Z,
								 int paramNum) {
		if (soln.atomicNumbers[atomNum1] != Z && soln.atomicNumbers[atomNum2] != Z) {
			return Ederivstandard(soln, atomNum1, atomNum2, densityDeriv, tau);
		}

		else if (paramNum != 1 && paramNum != 2 && paramNum != 5 && paramNum != 6) {
			return Ederivstandard(soln, atomNum1, atomNum2, densityDeriv, tau);
		}

		else if (paramNum == 1 || paramNum == 2) {//beta
			double e = Ederivstandard(soln, atomNum1, atomNum2, densityDeriv, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] orbsOfAtom = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			for (int i : orbsOfAtom[atomNum1]) {
				for (int k : orbsOfAtom[atomNum2]) {
					e += 2 * densityMatrix.get(i, k) *
							State.nom.Hbetapgd(orbitals[i], orbitals[k],
									ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
											soln.atomicNumbers[soln.atomOfOrb[k]], Z, orbitals[i].getL(),
											orbitals[k].getL(), paramNum - 1), tau);
				}
			}

			return e;
		}

		else {

			double e = Ederivstandard(soln, atomNum1, atomNum2, densityDeriv, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] orbsOfAtom = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : orbsOfAtom[atomNum1]) {
				for (int j : orbsOfAtom[atomNum1]) {
					e += densityMatrix.get(i, j) * atoms[atomNum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomNum1], soln.atomicNumbers[atomNum2], Z),
							paramNum - 5, tau);
				}
			}

			for (int k : orbsOfAtom[atomNum2]) {
				for (int l : orbsOfAtom[atomNum2]) {
					e -= densityMatrix.get(k, l) * atoms[atomNum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomNum2], soln.atomicNumbers[atomNum1], Z),
							paramNum - 5, tau);
				}
			}

			for (int i : orbsOfAtom[atomNum1]) {
				for (int k : orbsOfAtom[atomNum2]) {
					e += 2 * densityMatrix.get(i, k) * State.nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomNum1], soln.atomicNumbers[atomNum2], Z),
							paramNum - 5, tau);
				}
			}

			for (int i : orbsOfAtom[atomNum1]) {
				for (int j : orbsOfAtom[atomNum1]) {
					for (int k : orbsOfAtom[atomNum2]) {
						for (int l : orbsOfAtom[atomNum2]) {
							e += (densityMatrix.get(i, j) * densityMatrix.get(k, l) -
									densityMatrix.get(i, k) * 0.5 * densityMatrix.get(j, l))
									* State.nom.Gpgd(orbitals[i], orbitals[j], orbitals[k],
									orbitals[l], ParamDerivative.getNum(soln.atomicNumbers[atomNum1],
											soln.atomicNumbers[atomNum2], Z), paramNum - 5, tau);
						}
					}
				}
			}

			return e;
		}
	}

	private static double Ederiv(SolutionU soln, int atomNum1, int atomNum2, SimpleMatrix densityDerivAlpha,
								 SimpleMatrix densityDerivBeta, int tau, int Z, int paramNum) {

		if (soln.atomicNumbers[atomNum1] != Z && soln.atomicNumbers[atomNum2] != Z) {
			return Ederivstandard(soln, atomNum1, atomNum2, densityDerivAlpha, densityDerivBeta, tau);
		}

		else if (paramNum != 1 && paramNum != 2 && paramNum != 5 && paramNum != 6) {
			return Ederivstandard(soln, atomNum1, atomNum2, densityDerivAlpha, densityDerivBeta, tau);
		}

		else if (paramNum == 1 || paramNum == 2) {//beta
			double e = Ederivstandard(soln, atomNum1, atomNum2, densityDerivAlpha, densityDerivBeta, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] orbsOfAtom = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			for (int i : orbsOfAtom[atomNum1]) {
				for (int k : orbsOfAtom[atomNum2]) {
					e += 2 * densityMatrix.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z, orbitals[i].getL(),
									orbitals[k].getL(), paramNum - 1), tau);
				}
			}

			return e;
		}

		else {

			double e = Ederivstandard(soln, atomNum1, atomNum2, densityDerivAlpha, densityDerivBeta, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] orbsOfAtom = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : orbsOfAtom[atomNum1]) {
				for (int j : orbsOfAtom[atomNum1]) {
					e += densityMatrix.get(i, j) * atoms[atomNum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomNum1], soln.atomicNumbers[atomNum2], Z),
							paramNum - 5, tau);
				}
			}

			for (int k : orbsOfAtom[atomNum2]) {
				for (int l : orbsOfAtom[atomNum2]) {
					e -= densityMatrix.get(k, l) * atoms[atomNum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomNum2], soln.atomicNumbers[atomNum1], Z),
							paramNum - 5, tau);
				}
			}

			for (int i : orbsOfAtom[atomNum1]) {
				for (int k : orbsOfAtom[atomNum2]) {
					e += 2 * densityMatrix.get(i, k) * State.nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomNum1], soln.atomicNumbers[atomNum2], Z),
							paramNum - 5, tau);
				}
			}

			for (int i : orbsOfAtom[atomNum1]) {
				for (int j : orbsOfAtom[atomNum1]) {
					for (int k : orbsOfAtom[atomNum2]) {
						for (int l : orbsOfAtom[atomNum2]) {
							e += (densityMatrix.get(i, j) * densityMatrix.get(k, l) -
									soln.alphaDensity().get(i, k) * soln.alphaDensity().get(j, l) -
									soln.betaDensity().get(i, k) * soln.betaDensity().get(j, l))
									* State.nom.Gpgd(orbitals[i], orbitals[j], orbitals[k],
									orbitals[l], ParamDerivative.getNum(soln.atomicNumbers[atomNum1],
											soln.atomicNumbers[atomNum2], Z), paramNum - 5, tau);
						}
					}
				}
			}

			return e;
		}
	}
}
