package nddo.geometry;

import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import static nddo.State.nom;

public class GeometryDerivative {
	public static double gradient(NDDOAtom[] atoms, SolutionR soln, int atomnum, int tau) {
		SimpleMatrix densitymatrix = soln.densityMatrix();
		NDDOOrbital[] orbitals = soln.orbitals;
		int[][] index = soln.orbsOfAtom;
		int[] atomnumber = soln.atomOfOrb;
		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++)
							if (a != atomnum) sum += atoms[a].Vgd(orbitals[j], orbitals[k], tau);
					}
					else sum -= atoms[atomnum].Vgd(orbitals[j], orbitals[k], tau);
				}
				else if (atomnumber[j] == atomnum) sum += nom.Hgd(orbitals[j], orbitals[k], tau);
				else if (atomnumber[k] == atomnum) sum += nom.Hgd(orbitals[k], orbitals[j], tau);
				H.set(j, k, sum);
				H.set(k, j, sum);
			}
		SimpleMatrix G = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++)
							if (a != atomnum) for (int l : index[a])
								if (l > -1) for (int m : index[a])
									if (m > -1) sum += densitymatrix.get(l, m) *
											nom.Ggd(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
					}
					else for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnum])
							if (m > -1) sum += densitymatrix.get(l, m) *
									nom.Ggd(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
				}
				else if (atomnumber[j] == atomnum) {
					for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnumber[k]])
							if (m > -1) sum -= 0.5 * densitymatrix.get(l, m) *
									nom.Ggd(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
				}
				else if (atomnumber[k] == atomnum) for (int l : index[atomnum])
					if (l > -1) for (int m : index[atomnumber[j]])
						if (m > -1) sum -= 0.5 * densitymatrix.get(l, m) *
								nom.Ggd(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		SimpleMatrix F = H.copy().plus(G);
		double e = 0;
		for (int j = 0; j < orbitals.length; j++)
			for (int k = 0; k < orbitals.length; k++) e += 0.5 * densitymatrix.get(j, k) * (H.get(j, k) + F.get(j, k));
		for (int j = 0; j < atoms.length; j++) if (j != atomnum) e += atoms[atomnum].crfDeriv(atoms[j], tau);
		return e;
	}

	public static double gradientUnrestricted(NDDOAtom[] atoms, SolutionU soln, int atomnum, int tau) {
		SimpleMatrix alphadensity = soln.alphaDensity();
		SimpleMatrix betadensity = soln.betaDensity();
		NDDOOrbital[] orbitals = soln.orbitals;
		int[][] index = soln.orbsOfAtom;
		int[] atomnumber = soln.atomOfOrb;
		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++)
							if (a != atomnum) sum += atoms[a].Vgd(orbitals[j], orbitals[k], tau);
					}
					else sum -= atoms[atomnum].Vgd(orbitals[j], orbitals[k], tau);
				}
				else if (atomnumber[j] == atomnum) sum += nom.Hgd(orbitals[j], orbitals[k], tau);
				else if (atomnumber[k] == atomnum) sum += nom.Hgd(orbitals[k], orbitals[j], tau);
				H.set(j, k, sum);
				H.set(k, j, sum);
			}
		SimpleMatrix J = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) if (atomnumber[j] == atomnum) {
					for (int a = 0; a < atoms.length; a++)
						if (a != atomnum) for (int l : index[a])
							if (l > -1) for (int m : index[a])
								if (m > -1) sum += (alphadensity.get(l, m) + betadensity.get(l, m)) *
										nom.Ggd(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
				}
				else for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnum])
							if (m > -1) sum += (alphadensity.get(l, m) + betadensity.get(l, m)) *
									nom.Ggd(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
				J.set(j, k, sum);
				J.set(k, j, sum);
			}
		SimpleMatrix Ka = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] != atomnumber[k]) if (atomnumber[j] == atomnum) {
					for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnumber[k]])
							if (m > -1) sum -= alphadensity.get(l, m) *
									nom.Ggd(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
				}
				else if (atomnumber[k] == atomnum) for (int l : index[atomnum])
					if (l > -1) for (int m : index[atomnumber[j]])
						if (m > -1) sum -= alphadensity.get(l, m) *
								nom.Ggd(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
				Ka.set(j, k, sum);
				Ka.set(k, j, sum);
			}
		SimpleMatrix Kb = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] != atomnumber[k]) if (atomnumber[j] == atomnum) {
					for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnumber[k]])
							if (m > -1) sum -= betadensity.get(l, m) *
									nom.Ggd(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
				}
				else if (atomnumber[k] == atomnum) for (int l : index[atomnum])
					if (l > -1) for (int m : index[atomnumber[j]])
						if (m > -1) sum -= betadensity.get(l, m) *
								nom.Ggd(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
				Kb.set(j, k, sum);
				Kb.set(k, j, sum);
			}
		SimpleMatrix Fa = H.plus(J).plus(Ka);
		SimpleMatrix Fb = H.plus(J).plus(Kb);
		double e = 0;
		for (int j = 0; j < orbitals.length; j++)
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * alphadensity.get(j, k) * (H.get(j, k) + Fa.get(j, k));
				e += 0.5 * betadensity.get(j, k) * (H.get(j, k) + Fb.get(j, k));
			}
		for (int j = 0; j < atoms.length; j++) if (j != atomnum) e += atoms[atomnum].crfDeriv(atoms[j], tau);
		return e;
	}

	public static SimpleMatrix[][] gradientRoutine(SolutionR soln) {
		SimpleMatrix[] fockderivatives = new SimpleMatrix[soln.atoms.length * 3];
		SimpleMatrix grad = new SimpleMatrix(soln.atoms.length * 3, 1);
		int count = 0;
		for (int a = 0; a < soln.atoms.length; a++)
			for (int tau = 0; tau < 3; tau++) {
				SimpleMatrix[] matrices = staticderivs(soln.atoms, soln, a, tau);
				fockderivatives[count] = matrices[1];
				double sum = 0;
				for (int i = 0; i < matrices[1].numRows(); i++)
					for (int j = 0; j < matrices[1].numRows(); j++)
						sum += 0.5 * soln.densityMatrix().get(i, j) * (matrices[0].get(i, j) + matrices[1].get(i, j));
				for (int j = 0; j < soln.atoms.length; j++)
					if (j != a) sum += soln.atoms[a].crfDeriv(soln.atoms[j], tau);
				grad.set(count, 0, sum);
				count++;
			}
		return new SimpleMatrix[][]{new SimpleMatrix[]{grad}, fockderivatives};
	}

	public static SimpleMatrix[][] gradientRoutine(SolutionU soln) {
		SimpleMatrix[] fockderivativesalpha = new SimpleMatrix[soln.atoms.length * 3];
		SimpleMatrix[] fockderivativesbeta = new SimpleMatrix[soln.atoms.length * 3];
		SimpleMatrix grad = new SimpleMatrix(soln.atoms.length * 3, 1);
		int count = 0;
		for (int a = 0; a < soln.atoms.length; a++)
			for (int tau = 0; tau < 3; tau++) {
				SimpleMatrix[] matrices = staticderivs(soln.atoms, soln, a, tau);
				fockderivativesalpha[count] = matrices[1];
				fockderivativesbeta[count] = matrices[2];
				double sum = 0;
				for (int i = 0; i < matrices[1].numRows(); i++)
					for (int j = 0; j < matrices[1].numRows(); j++) {
						sum += 0.5 * soln.alphaDensity().get(i, j) * (matrices[0].get(i, j) + matrices[1].get(i, j));
						sum += 0.5 * soln.betaDensity().get(i, j) * (matrices[0].get(i, j) + matrices[2].get(i, j));
					}
				for (int j = 0; j < soln.atoms.length; j++)
					if (j != a) sum += soln.atoms[a].crfDeriv(soln.atoms[j], tau);
				grad.set(count, 0, sum);
				count++;
			}
		return new SimpleMatrix[][]{new SimpleMatrix[]{grad}, fockderivativesalpha, fockderivativesbeta};
	}

	public static double grad(SolutionR soln, int atomnum, int tau) {
		double e = 0;
		for (int a = 0; a < soln.atoms.length; a++)
			if (a != atomnum) {
				e += Ederiv(atomnum, a, soln.orbsOfAtom, soln.densityMatrix(), soln.atoms, soln.orbitals, tau);
				e += soln.atoms[atomnum].crfDeriv(soln.atoms[a], tau);
			}
		return e;
	}

	public static double grad(SolutionU soln, int atomnum, int tau) {
		double e = 0;
		for (int a = 0; a < soln.atoms.length; a++)
			if (a != atomnum) {
				e += Ederiv(atomnum, a, soln.orbsOfAtom, soln.alphaDensity(), soln.betaDensity(), soln.atoms,
						soln.orbitals, tau);
				e += soln.atoms[atomnum].crfDeriv(soln.atoms[a], tau);
			}
		return e;
	}

	private static double Ederiv(int atomnum1, int atomnum2, int[][] index, SimpleMatrix densityMatrix,
								 NDDOAtom[] atoms, NDDOOrbital[] orbitals, int tau) {
		double e = 0;
		for (int i : index[atomnum1])
			for (int j : index[atomnum1])
				if (i != -1 && j != -1)
					e += densityMatrix.get(i, j) * atoms[atomnum2].Vgd(orbitals[i], orbitals[j], tau);
		for (int k : index[atomnum2])
			for (int l : index[atomnum2])
				if (k != -1 && l != -1)
					e -= densityMatrix.get(k, l) * atoms[atomnum1].Vgd(orbitals[k], orbitals[l], tau);
		for (int i : index[atomnum1])
			for (int k : index[atomnum2])
				if (i != -1 && k != -1)
					e += 2 * densityMatrix.get(i, k) * nom.Hgd(orbitals[i], orbitals[k], tau);
		for (int i : index[atomnum1])
			for (int j : index[atomnum1])
				for (int k : index[atomnum2])
					for (int l : index[atomnum2])
						if (i != -1 && j != -1 && k != -1 && l != -1) e +=
								(densityMatrix.get(i, j) * densityMatrix.get(k, l) -
										densityMatrix.get(i, k) * 0.5 * densityMatrix.get(j, l)) *
										nom.Ggd(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau);
		return e;
	}

	private static double Ederiv(int atomnum1, int atomnum2, int[][] index, SimpleMatrix alphaDensity,
								 SimpleMatrix betaDensity, NDDOAtom[] atoms, NDDOOrbital[] orbitals, int tau) {
		double e = 0;
		for (int i : index[atomnum1])
			for (int j : index[atomnum1])
				if (i != -1 && j != -1) e += (alphaDensity.get(i, j) + betaDensity.get(i, j)) *
						atoms[atomnum2].Vgd(orbitals[i], orbitals[j], tau);
		for (int k : index[atomnum2])
			for (int l : index[atomnum2])
				if (k != -1 && l != -1) e -= (alphaDensity.get(k, l) + betaDensity.get(k, l)) *
						atoms[atomnum1].Vgd(orbitals[k], orbitals[l], tau);
		for (int i : index[atomnum1])
			for (int k : index[atomnum2])
				if (i != -1 && k != -1) e += 2 * (alphaDensity.get(i, k) + betaDensity.get(i, k)) *
						nom.Hgd(orbitals[i], orbitals[k], tau);
		for (int i : index[atomnum1])
			for (int j : index[atomnum1])
				for (int k : index[atomnum2])
					for (int l : index[atomnum2])
						if (i != -1 && j != -1 && k != -1 && l != -1) e +=
								((alphaDensity.get(i, j) + betaDensity.get(i, j)) *
										(alphaDensity.get(k, l) + betaDensity.get(k, l)) -
										alphaDensity.get(i, k) * alphaDensity.get(j, l) -
										betaDensity.get(i, k) * betaDensity.get(j, l)) *
										nom.Ggd(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau);
		return e;
	}

	public static SimpleMatrix[] staticderivs(NDDOAtom[] atoms, SolutionR soln, int atomnum, int tau) {
		SimpleMatrix densitymatrix = soln.densityMatrix();
		NDDOOrbital[] orbitals = soln.orbitals;
		int[][] index = soln.orbsOfAtom;
		int[] atomnumber = soln.atomOfOrb;
		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++)
							if (a != atomnum) sum += atoms[a].Vgd(orbitals[j], orbitals[k], tau);
					}
					else sum -= atoms[atomnum].Vgd(orbitals[j], orbitals[k], tau);
				}
				else if (atomnumber[j] == atomnum) sum += nom.Hgd(orbitals[j], orbitals[k], tau);
				else if (atomnumber[k] == atomnum) sum += nom.Hgd(orbitals[k], orbitals[j], tau);
				H.set(j, k, sum);
				H.set(k, j, sum);
			}
		SimpleMatrix G = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++)
							if (a != atomnum) for (int l : index[a])
								if (l > -1) for (int m : index[a])
									if (m > -1) sum += densitymatrix.get(l, m) *
											nom.Ggd(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
					}
					else for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnum])
							if (m > -1) sum += densitymatrix.get(l, m) *
									nom.Ggd(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
				}
				else if (atomnumber[j] == atomnum) {
					for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnumber[k]])
							if (m > -1) sum -= 0.5 * densitymatrix.get(l, m) *
									nom.Ggd(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
				}
				else if (atomnumber[k] == atomnum) for (int l : index[atomnum])
					if (l > -1) for (int m : index[atomnumber[j]])
						if (m > -1) sum -= 0.5 * densitymatrix.get(l, m) *
								nom.Ggd(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
				G.set(j, k, sum);
				G.set(k, j, sum);
			}
		SimpleMatrix F = H.copy().plus(G);
		return new SimpleMatrix[]{H, F};
	}

	public static SimpleMatrix[] staticderivs(NDDOAtom[] atoms, SolutionU soln, int atomnum, int tau) {
		SimpleMatrix alphadensity = soln.alphaDensity();
		SimpleMatrix betadensity = soln.betaDensity();
		NDDOOrbital[] orbitals = soln.orbitals;
		int[][] index = soln.orbsOfAtom;
		int[] atomnumber = soln.atomOfOrb;
		SimpleMatrix H = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++)
							if (a != atomnum) sum += atoms[a].Vgd(orbitals[j], orbitals[k], tau);
					}
					else sum -= atoms[atomnum].Vgd(orbitals[j], orbitals[k], tau);
				}
				else if (atomnumber[j] == atomnum) sum += nom.Hgd(orbitals[j], orbitals[k], tau);
				else if (atomnumber[k] == atomnum) sum += nom.Hgd(orbitals[k], orbitals[j], tau);
				H.set(j, k, sum);
				H.set(k, j, sum);
			}
		SimpleMatrix J = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) if (atomnumber[j] == atomnum) {
					for (int a = 0; a < atoms.length; a++)
						if (a != atomnum) for (int l : index[a])
							if (l > -1) for (int m : index[a])
								if (m > -1) sum += (alphadensity.get(l, m) + betadensity.get(l, m)) *
										nom.Ggd(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
				}
				else for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnum])
							if (m > -1) sum += (alphadensity.get(l, m) + betadensity.get(l, m)) *
									nom.Ggd(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
				J.set(j, k, sum);
				J.set(k, j, sum);
			}
		SimpleMatrix Ka = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] != atomnumber[k]) if (atomnumber[j] == atomnum) {
					for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnumber[k]])
							if (m > -1) sum -= alphadensity.get(l, m) *
									nom.Ggd(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
				}
				else if (atomnumber[k] == atomnum) for (int l : index[atomnum])
					if (l > -1) for (int m : index[atomnumber[j]])
						if (m > -1) sum -= alphadensity.get(l, m) *
								nom.Ggd(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
				Ka.set(j, k, sum);
				Ka.set(k, j, sum);
			}
		SimpleMatrix Kb = new SimpleMatrix(orbitals.length, orbitals.length);
		for (int j = 0; j < orbitals.length; j++)
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] != atomnumber[k]) if (atomnumber[j] == atomnum) {
					for (int l : index[atomnum])
						if (l > -1) for (int m : index[atomnumber[k]])
							if (m > -1) sum -= betadensity.get(l, m) *
									nom.Ggd(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
				}
				else if (atomnumber[k] == atomnum) for (int l : index[atomnum])
					if (l > -1) for (int m : index[atomnumber[j]])
						if (m > -1) sum -= betadensity.get(l, m) *
								nom.Ggd(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
				Kb.set(j, k, sum);
				Kb.set(k, j, sum);
			}
		SimpleMatrix Fa = H.plus(J).plus(Ka);
		SimpleMatrix Fb = H.plus(J).plus(Kb);
		return new SimpleMatrix[]{H, Fa, Fb};
	}

	public static SimpleMatrix densitymatrixderivfinite(NDDOAtom[] atoms, SolutionR soln, int atomnum, int tau) {
		SimpleMatrix orig = soln.densityMatrix();
		NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum, tau);
		SimpleMatrix perturbed = soln.withNewAtoms(newatoms).densityMatrix();
		return perturbed.minus(orig).scale(1E7);
	}

	public static SimpleMatrix[] densitymatrixderivfinite(NDDOAtom[] atoms, SolutionU soln, int atomnum, int tau) {
		SimpleMatrix aorig = soln.alphaDensity();
		SimpleMatrix borig = soln.betaDensity();
		NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum, tau);
		SolutionU newsoln = (SolutionU) soln.withNewAtoms(newatoms);
		SimpleMatrix aperturbed = newsoln.alphaDensity();
		SimpleMatrix bperturbed = newsoln.betaDensity();
		return new SimpleMatrix[]{aperturbed.minus(aorig).scale(1E7), bperturbed.minus(borig).scale(1E7)};
	}
}
