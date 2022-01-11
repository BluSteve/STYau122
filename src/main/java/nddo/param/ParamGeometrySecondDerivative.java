package nddo.param;

import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.solution.Solution;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;

import java.util.stream.IntStream;

import static nddo.State.nom;

public class ParamGeometrySecondDerivative {
	public static double gradDeriv2Alpha(Solution soln, int atomNum, int tau, int Z) {
		double e = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			if (a != atomNum) {
				e += soln.atoms[atomNum].crfalphap2gd(soln.atoms[a],
						ParamDerivative.getNum(soln.atoms[atomNum].getAtomProperties().getZ(),
								soln.atoms[a].getAtomProperties().getZ(), Z), tau);
			}
		}

		return e;
	}

	public static double gradderiv2(SolutionR soln, int atomnum, int tau, int Z1, int paramnum1, int Z2, int paramnum2,
									SimpleMatrix densityderiva, SimpleMatrix densityderivb,
									SimpleMatrix densityderiv2) {
		return IntStream.range(0, soln.atoms.length).parallel().mapToDouble(
				a -> {
					if (a != atomnum) {
						return Ederiv2(soln, atomnum, a, densityderiva, densityderivb, densityderiv2, tau, Z1,
								paramnum1, Z2, paramnum2);
					}

					return 0;
				}).sum();
	}

	public static double gradderiv2(SolutionU soln, int atomnum, int tau, int Z1, int paramnum1, int Z2, int paramnum2,
									SimpleMatrix densityderivaalpha, SimpleMatrix densityderivabeta,
									SimpleMatrix densityderivbalpha, SimpleMatrix densityderivbbeta,
									SimpleMatrix densityderiv2alpha, SimpleMatrix densityderiv2beta) {

		return IntStream.range(0, soln.atoms.length).parallel().mapToDouble(
				a -> {
					if (a != atomnum) {
						return Ederiv2(soln, atomnum, a, densityderivaalpha, densityderivabeta, densityderivbalpha,
								densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau, Z1, paramnum1, Z2,
								paramnum2);
					}

					return 0;
				}).sum();
	}

	private static double Ederiv2(SolutionR soln, int atomnum1, int atomnum2, SimpleMatrix densityderiva,
								  SimpleMatrix densityderivb, SimpleMatrix densityderiv2, int tau, int Z1,
								  int paramnum1, int Z2, int paramnum2) {

		if (soln.atomicNumbers[atomnum1] != Z1 && soln.atomicNumbers[atomnum2] != Z1 &&
				soln.atomicNumbers[atomnum1] != Z2 && soln.atomicNumbers[atomnum2] != Z2) {
			return Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);
		}

		else if (paramnum1 != 1 && paramnum1 != 2 && paramnum1 != 5 && paramnum1 != 6 && paramnum2 != 1 &&
				paramnum2 != 2 && paramnum2 != 5 && paramnum2 != 6) {
			return Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);
		}

		else if ((paramnum1 == 1 || paramnum1 == 2) && paramnum2 != 1 && paramnum2 != 2 && paramnum2 != 5 &&
				paramnum2 != 6) {//beta-irrelevant

			double e = Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderivb.get(i, k) * nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z1, orbitals[i].getL(), orbitals[k].getL(),
									paramnum1 - 1), tau);
				}
			}

			return e;
		}

		else if ((paramnum2 == 1 || paramnum2 == 2) && paramnum1 != 1 && paramnum1 != 2 && paramnum1 != 5 &&
				paramnum1 != 6) {//irrelevant-beta

			double e = Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z2, orbitals[i].getL(), orbitals[k].getL(),
									paramnum2 - 1), tau);
				}
			}

			return e;
		}

		else if ((paramnum2 == 1 || paramnum2 == 2) && (paramnum1 == 1 || paramnum1 == 2)) {//beta-beta

			double e = Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z2, orbitals[i].getL(), orbitals[k].getL(),
									paramnum2 - 1), tau);

					e += 2 * densityderivb.get(i, k) * nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z1, orbitals[i].getL(), orbitals[k].getL(),
									paramnum1 - 1), tau);
				}
			}

			return e;
		}

		else if ((paramnum1 == 5 || paramnum1 == 6) && paramnum2 != 1 && paramnum2 != 2 && paramnum2 != 5 &&
				paramnum2 != 6) {//zeta-irrelevant

			double e = Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderivb.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderivb.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderivb.get(i, k) * nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderivb.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderivb.get(k, l) -
									densityderivb.get(i, k) * 0.5 * densityMatrix.get(j, l) -
									densityMatrix.get(i, k) * 0.5 * densityderivb.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z1), paramnum1 - 5, tau);
						}
					}
				}
			}

			return e;
		}

		else if ((paramnum2 == 5 || paramnum2 == 6) && paramnum1 != 1 && paramnum1 != 2 && paramnum1 != 5 &&
				paramnum1 != 6) {//irrelevant-zeta

			double e = Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderiva.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderiva.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderiva.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderiva.get(k, l) -
									densityderiva.get(i, k) * 0.5 * densityMatrix.get(j, l) -
									densityMatrix.get(i, k) * 0.5 * densityderiva.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z2), paramnum2 - 5, tau);
						}
					}
				}
			}

			return e;
		}

		else if ((paramnum2 == 5 || paramnum2 == 6) && (paramnum1 == 1 || paramnum1 == 2)) {//beta-zeta

			double e = Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderiva.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderiva.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * State.nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);

					e += 2 * densityderivb.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z1, orbitals[i].getL(), orbitals[k].getL(),
									paramnum1 - 1), tau);

					e += 2 * densityMatrix.get(i, k) * State.nom.Hbetazetap2gd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z1, orbitals[i].getL(), orbitals[k].getL(),
									paramnum1 - 1),
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderiva.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderiva.get(k, l) -
									densityderiva.get(i, k) * 0.5 * densityMatrix.get(j, l) -
									densityMatrix.get(i, k) * 0.5 * densityderiva.get(j, l)) *
									State.nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z2), paramnum2 - 5, tau);
						}
					}
				}
			}

			return e;
		}

		else if ((paramnum1 == 5 || paramnum1 == 6) && (paramnum2 == 1 || paramnum2 == 2)) {//zeta-beta

			double e = Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderivb.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderivb.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderivb.get(i, k) * State.nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);

					e += 2 * densityderiva.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z2, orbitals[i].getL(), orbitals[k].getL(),
									paramnum2 - 1), tau);

					e += 2 * densityMatrix.get(i, k) * State.nom.Hbetazetap2gd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z2, orbitals[i].getL(), orbitals[k].getL(),
									paramnum2 - 1),
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderivb.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderivb.get(k, l) -
									densityderivb.get(i, k) * 0.5 * densityMatrix.get(j, l) -
									densityMatrix.get(i, k) * 0.5 * densityderivb.get(j, l)) *
									State.nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z1), paramnum1 - 5, tau);
						}
					}
				}
			}

			return e;
		}

		else if ((paramnum1 == 5 || paramnum1 == 6) && (paramnum2 == 5 || paramnum2 == 6)) {//zeta-zeta

			double e = Ederiv2standard(soln, atomnum1, atomnum2, densityderiva, densityderivb, densityderiv2, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderiva.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
					e += densityderivb.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
					e += densityMatrix.get(i, j) * atoms[atomnum2].Vp2gd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5,
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderiva.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z2),
							paramnum2 - 5, tau);
					e -= densityderivb.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z1),
							paramnum1 - 5, tau);
					e -= densityMatrix.get(k, l) * atoms[atomnum1].Vp2gd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z1),
							paramnum1 - 5,
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * State.nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
					e += 2 * densityderivb.get(i, k) * State.nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
					e += 2 * densityMatrix.get(i, k) * State.nom.Hzetazetap2gd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5,
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderiva.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderiva.get(k, l) -
									densityderiva.get(i, k) * 0.5 * densityMatrix.get(j, l) -
									densityMatrix.get(i, k) * 0.5 * densityderiva.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z2), paramnum2 - 5, tau);

							e += (densityderivb.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderivb.get(k, l) -
									densityderivb.get(i, k) * 0.5 * densityMatrix.get(j, l) -
									densityMatrix.get(i, k) * 0.5 * densityderivb.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z1), paramnum1 - 5, tau);

							e += (densityMatrix.get(i, j) * densityMatrix.get(k, l) -
									densityMatrix.get(i, k) * 0.5 * densityMatrix.get(j, l)) *
									State.nom.Gp2gd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z1), paramnum1 - 5,
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z2), paramnum2 - 5, tau);
						}
					}
				}
			}

			return e;

		}

		else {
			System.err.println("not done yet");
			return 0;
		}


	}

	private static double Ederiv2standard(SolutionR soln, int atomnum1, int atomnum2, SimpleMatrix densityderiva,
										  SimpleMatrix densityderivb, SimpleMatrix densityderiv2, int tau) {

		int[][] index = soln.orbsOfAtom;

		SimpleMatrix densityMatrix = soln.densityMatrix();

		NDDOAtom[] atoms = soln.atoms;

		NDDOOrbital[] orbitals = soln.orbitals;

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				e += densityderiv2.get(i, j) * atoms[atomnum2].Vgd(orbitals[i], orbitals[j], tau);
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderiv2.get(k, l) * atoms[atomnum1].Vgd(orbitals[k], orbitals[l], tau);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				e += 2 * densityderiv2.get(i, k) * nom.Hgd(orbitals[i], orbitals[k], tau);
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						e += (densityderiv2.get(i, j) * densityMatrix.get(k, l) +
								densityMatrix.get(i, j) * densityderiv2.get(k, l) -
								densityderiv2.get(i, k) * 0.5 * densityMatrix.get(j, l) -
								densityMatrix.get(i, k) * 0.5 * densityderiv2.get(j, l) +
								densityderiva.get(i, j) * densityderivb.get(k, l) +
								densityderivb.get(i, j) * densityderiva.get(k, l) -
								densityderiva.get(i, k) * 0.5 * densityderivb.get(j, l) -
								densityderivb.get(i, k) * 0.5 * densityderiva.get(j, l)) *
								nom.Ggd(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau);
					}
				}
			}
		}

		return e;
	}

	private static double Ederiv2(SolutionU soln, int atomnum1, int atomnum2, SimpleMatrix densityderivaalpha,
								  SimpleMatrix densityderivabeta, SimpleMatrix densityderivbalpha,
								  SimpleMatrix densityderivbbeta, SimpleMatrix densityderiv2alpha,
								  SimpleMatrix densityderiv2beta, int tau, int Z1, int paramnum1, int Z2,
								  int paramnum2) {

		if (soln.atomicNumbers[atomnum1] != Z1 && soln.atomicNumbers[atomnum2] != Z1 &&
				soln.atomicNumbers[atomnum1] != Z2 && soln.atomicNumbers[atomnum2] != Z2) {
			return Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta, densityderivbalpha,
					densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);
		}

		else if (paramnum1 != 1 && paramnum1 != 2 && paramnum1 != 5 && paramnum1 != 6 && paramnum2 != 1 &&
				paramnum2 != 2 && paramnum2 != 5 && paramnum2 != 6) {
			return Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta, densityderivbalpha,
					densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);
		}

		else if ((paramnum1 == 1 || paramnum1 == 2) && paramnum2 != 1 && paramnum2 != 2 && paramnum2 != 5 &&
				paramnum2 != 6) {//beta-irrelevant

			double e =
					Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta,
							densityderivbalpha,
							densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);

			SimpleMatrix densityderivb = densityderivbalpha.plus(densityderivbbeta);

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderivb.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z1, orbitals[i].getL(), orbitals[k].getL(),
									paramnum1 - 1), tau);
				}
			}

			return e;
		}

		else if ((paramnum2 == 1 || paramnum2 == 2) && paramnum1 != 1 && paramnum1 != 2 && paramnum1 != 5 &&
				paramnum1 != 6) {//irrelevant-beta

			double e =
					Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta,
							densityderivbalpha,
							densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);

			SimpleMatrix densityderiva = densityderivaalpha.plus(densityderivabeta);

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z2, orbitals[i].getL(), orbitals[k].getL(),
									paramnum2 - 1), tau);
				}
			}

			return e;
		}

		else if ((paramnum2 == 1 || paramnum2 == 2) && (paramnum1 == 1 || paramnum1 == 2)) {//beta-beta

			double e =
					Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta,
							densityderivbalpha,
							densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			SimpleMatrix densityderiva = densityderivaalpha.plus(densityderivabeta);

			SimpleMatrix densityderivb = densityderivbalpha.plus(densityderivbbeta);


			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z2, orbitals[i].getL(), orbitals[k].getL(),
									paramnum2 - 1), tau);

					e += 2 * densityderivb.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z1, orbitals[i].getL(), orbitals[k].getL(),
									paramnum1 - 1), tau);
				}
			}

			return e;
		}

		else if ((paramnum1 == 5 || paramnum1 == 6) && paramnum2 != 1 && paramnum2 != 2 && paramnum2 != 5 &&
				paramnum2 != 6) {//zeta-irrelevant

			double e =
					Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta,
							densityderivbalpha,
							densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			SimpleMatrix densityderivb = densityderivbalpha.plus(densityderivbbeta);

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderivb.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderivb.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderivb.get(i, k) * State.nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderivb.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderivb.get(k, l) -
									densityderivbalpha.get(i, k) * soln.alphaDensity().get(j, l) -
									densityderivbbeta.get(i, k) * soln.betaDensity().get(j, l) -
									soln.alphaDensity().get(i, k) * densityderivbalpha.get(j, l) -
									soln.betaDensity().get(i, k) * densityderivbbeta.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z1), paramnum1 - 5, tau);
						}
					}
				}
			}

			return e;
		}

		else if ((paramnum2 == 5 || paramnum2 == 6) && paramnum1 != 1 && paramnum1 != 2 && paramnum1 != 5 &&
				paramnum1 != 6) {//irrelevant-zeta

			double e =
					Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta,
							densityderivbalpha,
							densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			SimpleMatrix densityderiva = densityderivaalpha.plus(densityderivabeta);


			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderiva.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderiva.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderiva.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderiva.get(k, l) -
									densityderivaalpha.get(i, k) * soln.alphaDensity().get(j, l) -
									densityderivabeta.get(i, k) * soln.betaDensity().get(j, l) -
									soln.alphaDensity().get(i, k) * densityderivaalpha.get(j, l) -
									soln.betaDensity().get(i, k) * densityderivabeta.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z2), paramnum2 - 5, tau);
						}
					}
				}
			}

			return e;
		}

		else if ((paramnum2 == 5 || paramnum2 == 6) && (paramnum1 == 1 || paramnum1 == 2)) {//beta-zeta

			double e =
					Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta,
							densityderivbalpha,
							densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			SimpleMatrix densityderiva = densityderivaalpha.plus(densityderivabeta);

			SimpleMatrix densityderivb = densityderivbalpha.plus(densityderivbbeta);

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderiva.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderiva.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);

					e += 2 * densityderivb.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z1, orbitals[i].getL(), orbitals[k].getL(),
									paramnum1 - 1), tau);

					e += 2 * densityMatrix.get(i, k) * State.nom.Hbetazetap2gd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z1, orbitals[i].getL(), orbitals[k].getL(),
									paramnum1 - 1),
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderiva.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderiva.get(k, l) -
									densityderivaalpha.get(i, k) * soln.alphaDensity().get(j, l) -
									densityderivabeta.get(i, k) * soln.betaDensity().get(j, l) -
									soln.alphaDensity().get(i, k) * densityderivaalpha.get(j, l) -
									soln.betaDensity().get(i, k) * densityderivabeta.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z2), paramnum2 - 5, tau);
						}
					}
				}
			}

			return e;
		}

		else if ((paramnum1 == 5 || paramnum1 == 6) && (paramnum2 == 1 || paramnum2 == 2)) {//zeta-beta

			double e =
					Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta,
							densityderivbalpha,
							densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			SimpleMatrix densityderiva = densityderivaalpha.plus(densityderivabeta);

			SimpleMatrix densityderivb = densityderivbalpha.plus(densityderivbbeta);

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderivb.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderivb.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderivb.get(i, k) * nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);

					e += 2 * densityderiva.get(i, k) * State.nom.Hbetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z2, orbitals[i].getL(), orbitals[k].getL(),
									paramnum2 - 1), tau);

					e += 2 * densityMatrix.get(i, k) * State.nom.Hbetazetap2gd(orbitals[i], orbitals[k],
							ParamDerivative.getNumBeta(soln.atomicNumbers[soln.atomOfOrb[i]],
									soln.atomicNumbers[soln.atomOfOrb[k]], Z2, orbitals[i].getL(), orbitals[k].getL(),
									paramnum2 - 1),
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderivb.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderivb.get(k, l) -
									densityderivbalpha.get(i, k) * soln.alphaDensity().get(j, l) -
									densityderivbbeta.get(i, k) * soln.betaDensity().get(j, l) -
									soln.alphaDensity().get(i, k) * densityderivbalpha.get(j, l) -
									soln.betaDensity().get(i, k) * densityderivbbeta.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z1), paramnum1 - 5, tau);
						}
					}
				}
			}

			return e;
		}

		else if ((paramnum1 == 5 || paramnum1 == 6) && (paramnum2 == 5 || paramnum2 == 6)) {//zeta-zeta

			double e =
					Ederiv2standard(soln, atomnum1, atomnum2, densityderivaalpha, densityderivabeta,
							densityderivbalpha,
							densityderivbbeta, densityderiv2alpha, densityderiv2beta, tau);

			SimpleMatrix densityMatrix = soln.densityMatrix();

			SimpleMatrix densityderiva = densityderivaalpha.plus(densityderivabeta);

			SimpleMatrix densityderivb = densityderivbalpha.plus(densityderivbbeta);

			int[][] index = soln.orbsOfAtom;

			NDDOOrbital[] orbitals = soln.orbitals;

			NDDOAtom[] atoms = soln.atoms;

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					e += densityderiva.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
					e += densityderivb.get(i, j) * atoms[atomnum2].Vpgd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
					e += densityMatrix.get(i, j) * atoms[atomnum2].Vp2gd(orbitals[i], orbitals[j],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5,
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int k : index[atomnum2]) {
				for (int l : index[atomnum2]) {
					e -= densityderiva.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z2),
							paramnum2 - 5, tau);
					e -= densityderivb.get(k, l) * atoms[atomnum1].Vpgd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z1),
							paramnum1 - 5, tau);
					e -= densityMatrix.get(k, l) * atoms[atomnum1].Vp2gd(orbitals[k], orbitals[l],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z1),
							paramnum1 - 5,
							ParamDerivative.getNum(soln.atomicNumbers[atomnum2], soln.atomicNumbers[atomnum1], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					e += 2 * densityderiva.get(i, k) * nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
					e += 2 * densityderivb.get(i, k) * nom.Hzetapgd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5, tau);
					e += 2 * densityMatrix.get(i, k) * nom.Hzetazetap2gd(orbitals[i], orbitals[k],
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z1),
							paramnum1 - 5,
							ParamDerivative.getNum(soln.atomicNumbers[atomnum1], soln.atomicNumbers[atomnum2], Z2),
							paramnum2 - 5, tau);
				}
			}

			for (int i : index[atomnum1]) {
				for (int j : index[atomnum1]) {
					for (int k : index[atomnum2]) {
						for (int l : index[atomnum2]) {
							e += (densityderiva.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderiva.get(k, l) -
									densityderivaalpha.get(i, k) * soln.alphaDensity().get(j, l) -
									densityderivabeta.get(i, k) * soln.betaDensity().get(j, l) -
									soln.alphaDensity().get(i, k) * densityderivaalpha.get(j, l) -
									soln.betaDensity().get(i, k) * densityderivabeta.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z2), paramnum2 - 5, tau);

							e += (densityderivb.get(i, j) * densityMatrix.get(k, l) +
									densityMatrix.get(i, j) * densityderivb.get(k, l) -
									densityderivbalpha.get(i, k) * soln.alphaDensity().get(j, l) -
									densityderivbbeta.get(i, k) * soln.betaDensity().get(j, l) -
									soln.alphaDensity().get(i, k) * densityderivbalpha.get(j, l) -
									soln.betaDensity().get(i, k) * densityderivbbeta.get(j, l)) *
									nom.Gpgd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z1), paramnum1 - 5, tau);

							e += (densityMatrix.get(i, j) * densityMatrix.get(k, l) -
									soln.alphaDensity().get(i, k) * soln.alphaDensity().get(j, l) -
									soln.betaDensity().get(i, k) * soln.betaDensity().get(j, l)) *
									State.nom.Gp2gd(orbitals[i], orbitals[j], orbitals[k], orbitals[l],
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z1), paramnum1 - 5,
											ParamDerivative.getNum(soln.atomicNumbers[atomnum1],
													soln.atomicNumbers[atomnum2], Z2), paramnum2 - 5, tau);
						}
					}
				}
			}

			return e;

		}

		else {
			System.err.println("not done yet");
			return 0;
		}


	}

	private static double Ederiv2standard(SolutionU soln, int atomnum1, int atomnum2, SimpleMatrix densityderivaalpha,
										  SimpleMatrix densityderivabeta, SimpleMatrix densityderivbalpha,
										  SimpleMatrix densityderivbbeta, SimpleMatrix densityderiv2alpha,
										  SimpleMatrix densityderiv2beta, int tau) {

		int[][] index = soln.orbsOfAtom;

		SimpleMatrix densityMatrix = soln.densityMatrix();

		NDDOAtom[] atoms = soln.atoms;

		NDDOOrbital[] orbitals = soln.orbitals;

		SimpleMatrix densityderiva = densityderivaalpha.plus(densityderivabeta);

		SimpleMatrix densityderivb = densityderivbalpha.plus(densityderivbbeta);

		SimpleMatrix densityderiv2 = densityderiv2alpha.plus(densityderiv2beta);

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				e += densityderiv2.get(i, j) * atoms[atomnum2].Vgd(orbitals[i], orbitals[j], tau);
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				e -= densityderiv2.get(k, l) * atoms[atomnum1].Vgd(orbitals[k], orbitals[l], tau);
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				e += 2 * densityderiv2.get(i, k) * nom.Hgd(orbitals[i], orbitals[k], tau);
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						e += (densityderiv2.get(i, j) * densityMatrix.get(k, l) +
								densityMatrix.get(i, j) * densityderiv2.get(k, l) -
								densityderiv2alpha.get(i, k) * soln.alphaDensity().get(j, l) -
								densityderiv2beta.get(i, k) * soln.betaDensity().get(j, l) -
								soln.alphaDensity().get(i, k) * densityderiv2alpha.get(j, l) -
								soln.betaDensity().get(i, k) * densityderiv2beta.get(j, l) +
								densityderiva.get(i, j) * densityderivb.get(k, l) +
								densityderivb.get(i, j) * densityderiva.get(k, l) -
								densityderivaalpha.get(i, k) * densityderivbalpha.get(j, l) -
								densityderivabeta.get(i, k) * densityderivbbeta.get(j, l) -
								densityderivbalpha.get(i, k) * densityderivaalpha.get(j, l) -
								densityderivbbeta.get(i, k) * densityderivabeta.get(j, l)) *
								nom.Ggd(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau);
					}
				}
			}
		}

		return e;
	}
}
