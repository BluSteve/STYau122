package nddo.geometry;

import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.math.PopleThiel;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Batcher;

import static nddo.State.nom;

public class GeometrySecondDerivative {
	public static SimpleMatrix hessianRoutine(SolutionR soln, SimpleMatrix[] fockderivstatic) {
		SimpleMatrix[] densityDerivs;
		if (soln.rm.nonvAlpha == 0) {
			SimpleMatrix sm = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
			densityDerivs = new SimpleMatrix[fockderivstatic.length];
			for (int i = 0; i < fockderivstatic.length; i++) {
				densityDerivs[i] = sm;
			}
		}
		else densityDerivs = Batcher.apply(fockderivstatic, State.config.poplethiel_batch_size, subset -> densityDeriv(soln, subset));

		SimpleMatrix hessian = new SimpleMatrix(densityDerivs.length, densityDerivs.length);

		int[][] indices = new int[hessian.numRows() * (hessian.numRows() + 1) / 2][];
		int count = 0;
		for (int i = 0; i < hessian.numRows(); i++) {
			for (int j = i; j < hessian.numCols(); j++) {
				indices[count] = new int[]{i, j};
				count++;
			}
		}

		Batcher.consume(indices, 1, subset -> {
			for (int[] ints : subset) {
				int i = ints[0];
				int j = ints[1];
				double E = 0;
				int atomnum1 = i / 3;
				int atomnum2 = j / 3;
				int tau1 = i - 3 * atomnum1;
				int tau2 = j - 3 * atomnum2;

				if (atomnum1 == atomnum2) {
					for (int a = 0; a < soln.atoms.length; a++) {
						if (a != atomnum1) {
							E += Ederiv2(atomnum1, a, soln.orbsOfAtom, soln.densityMatrix(), soln.atoms, soln.orbitals,
									tau1, tau2);
							E += soln.atoms[atomnum1].crfg2d(soln.atoms[a], tau1, tau2);
						}
					}
				}
				else E = -Ederiv2(atomnum1, atomnum2, soln.orbsOfAtom, soln.densityMatrix(), soln.atoms, soln.orbitals,
						tau1, tau2) - soln.atoms[atomnum1].crfg2d(soln.atoms[atomnum2], tau1, tau2);

				for (int I = 0; I < soln.orbitals.length; I++) {
					for (int J = 0; J < soln.orbitals.length; J++) {
						E += fockderivstatic[i].get(I, J) * densityDerivs[j].get(I, J);
					}
				}

				hessian.set(i, j, E);
				hessian.set(j, i, E);
			}
		});

		return hessian;
	}

	public static SimpleMatrix hessianRoutine(SolutionU soln, SimpleMatrix[] fockderivstaticalpha,
											  SimpleMatrix[] fockderivstaticbeta) {
		SimpleMatrix[][] sms;
		if (soln.rm.nonvAlpha + soln.rm.nonvBeta == 0) {
			SimpleMatrix sm = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
			sms = new SimpleMatrix[2][fockderivstaticalpha.length];
			for (int i = 0; i < fockderivstaticalpha.length; i++) {
				sms[0][i] = sm;
				sms[1][i] = sm;
			}
		}
		else sms = Batcher.apply(fockderivstaticalpha, fockderivstaticbeta, SimpleMatrix[][].class, State.config.poplethiel_batch_size,
				(subseta, subsetb) -> densityDeriv(soln, subseta, subsetb));

		SimpleMatrix[] densityderivsalpha = sms[0];
		SimpleMatrix[] densityderivsbeta = sms[1];

		SimpleMatrix hessian = new SimpleMatrix(densityderivsalpha.length, densityderivsalpha.length);

		int[][] indices = new int[hessian.numRows() * (hessian.numRows() + 1) / 2][];
		int count = 0;
		for (int i = 0; i < hessian.numRows(); i++) {
			for (int j = i; j < hessian.numCols(); j++) {
				indices[count] = new int[]{i, j};
				count++;
			}
		}

		Batcher.consume(indices, 1, subset -> {
			for (int[] ints : subset) {
				int i = ints[0];
				int j = ints[1];
				double E = 0;
				int atomnum1 = i / 3;
				int atomnum2 = j / 3;
				int tau1 = i - 3 * atomnum1;
				int tau2 = j - 3 * atomnum2;

				if (atomnum1 == atomnum2) {
					for (int a = 0; a < soln.atoms.length; a++) {
						if (a != atomnum1) {
							E += Ederiv2(atomnum1, a, soln.orbsOfAtom, soln.alphaDensity(), soln.betaDensity(),
									soln.atoms, soln.orbitals, tau1, tau2);
							E += soln.atoms[atomnum1].crfg2d(soln.atoms[a], tau1, tau2);
						}
					}
				}
				else {
					E = -Ederiv2(atomnum1, atomnum2, soln.orbsOfAtom, soln.alphaDensity(), soln.betaDensity(),
							soln.atoms, soln.orbitals, tau1, tau2) -
							soln.atoms[atomnum1].crfg2d(soln.atoms[atomnum2], tau1, tau2);
				}

				for (int I = 0; I < soln.orbitals.length; I++) {
					for (int J = 0; J < soln.orbitals.length; J++) {
						E += fockderivstaticalpha[i].get(I, J) * densityderivsalpha[j].get(I, J);
						E += fockderivstaticbeta[i].get(I, J) * densityderivsbeta[j].get(I, J);
					}
				}

				hessian.set(i, j, E);
				hessian.set(j, i, E);
			}
		});

		return hessian;
	}

	private static SimpleMatrix[] densityDeriv(SolutionR soln, SimpleMatrix[] fockderivstatic) {
		SimpleMatrix[] xarray = PopleThiel.pople(soln, PopleThiel.toMO(soln.CtOcc, soln.CVirt, fockderivstatic));

		SimpleMatrix[] densityDerivs = new SimpleMatrix[fockderivstatic.length];

		for (int i = 0; i < xarray.length; i++) {
			densityDerivs[i] = PopleThiel.densityDeriv(soln, xarray[i]);
		}

		return densityDerivs;
	}

	private static SimpleMatrix[][] densityDeriv(SolutionU soln, SimpleMatrix[] fockderivstaticalpha, SimpleMatrix[] fockderivstaticbeta) {
		SimpleMatrix[] xarray = PopleThiel.thiel(soln,
				PopleThiel.toMO(soln.CtaOcc, soln.CaVirt, fockderivstaticalpha),
				PopleThiel.toMO(soln.CtbOcc, soln.CbVirt, fockderivstaticbeta));

		SimpleMatrix[] densityDerivsAlpha = new SimpleMatrix[fockderivstaticalpha.length];
		SimpleMatrix[] densityDerivsBeta = new SimpleMatrix[fockderivstaticbeta.length];

		for (int i = 0; i < xarray.length; i++) {
			SimpleMatrix[] sms = PopleThiel.densityDeriv(soln, xarray[i]);

			densityDerivsAlpha[i] = sms[0];
			densityDerivsBeta[i] = sms[1];
		}

		return new SimpleMatrix[][]{densityDerivsAlpha, densityDerivsBeta};
	}

	private static double Ederiv2(int atomnum1, int atomnum2, int[][] index, SimpleMatrix densityMatrix,
								  NDDOAtom[] atoms, NDDOOrbital[] orbitals, int tau1, int tau2) {
		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				e += densityMatrix.get(i, j) * atoms[atomnum2].Vg2d(orbitals[i], orbitals[j], tau1, tau2);
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				e += densityMatrix.get(k, l) * atoms[atomnum1].Vg2d(orbitals[k], orbitals[l], tau1, tau2);
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				e += 2 * densityMatrix.get(i, k) * nom.Hg2d(orbitals[i], orbitals[k], tau1, tau2);
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						e += (densityMatrix.get(i, j) * densityMatrix.get(k, l) -
								densityMatrix.get(i, k) * 0.5 * densityMatrix.get(j, l)) *
								nom.Gg2d(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau1, tau2);
					}
				}
			}
		}

		return e;
	}

	private static double Ederiv2(int atomnum1, int atomnum2, int[][] index, SimpleMatrix alphaDensity,
								  SimpleMatrix betaDensity, NDDOAtom[] atoms, NDDOOrbital[] orbitals, int tau1,
								  int tau2) {
		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				e += (alphaDensity.get(i, j) + betaDensity.get(i, j)) *
						atoms[atomnum2].Vg2d(orbitals[i], orbitals[j], tau1, tau2);
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				e += (alphaDensity.get(k, l) + betaDensity.get(k, l)) *
						atoms[atomnum1].Vg2d(orbitals[k], orbitals[l], tau1, tau2);
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				e += 2 * (alphaDensity.get(i, k) + betaDensity.get(i, k)) *
						nom.Hg2d(orbitals[i], orbitals[k], tau1, tau2);
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						e += ((alphaDensity.get(i, j) + betaDensity.get(i, j)) *
								(alphaDensity.get(k, l) + betaDensity.get(k, l)) -
								alphaDensity.get(i, k) * alphaDensity.get(j, l) -
								betaDensity.get(i, k) * betaDensity.get(j, l)) *
								nom.Gg2d(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau1, tau2);
					}
				}
			}
		}

		return e;
	}
}
