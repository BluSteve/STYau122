package nddo.geometry;

import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.State;
import nddo.math.PopleThiel;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.simple.SimpleMatrix;
import tools.Batcher;

import java.util.stream.IntStream;

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
		else densityDerivs = Batcher.apply(fockderivstatic, State.config.poplethiel_batch_size,
				subset -> {
					SimpleMatrix[] sms = PopleThiel.pt(soln, PopleThiel.toMO(soln.CtOcc, soln.CVirt, subset));
					SimpleMatrix[] results = new SimpleMatrix[sms.length];

					for (int j = 0; j < sms.length; j++) {
						results[j] = PopleThiel.densityDeriv(soln, sms[j]);
					}

					return results;
				});

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
					E += IntStream.range(0, soln.atoms.length).parallel().mapToDouble(a -> {
						double E2 = 0;

						if (a != atomnum1) {
							E2 += Ederiv2(atomnum1, a, soln.orbsOfAtom, soln.densityMatrix(),
									soln.atoms, soln.orbitals, tau1, tau2);
							E2 += soln.atoms[atomnum1].crfg2d(soln.atoms[a], tau1, tau2);
						}

						return E2;
					}).sum();
				}
				else E = -Ederiv2(atomnum1, atomnum2, soln.orbsOfAtom, soln.densityMatrix(),
						soln.atoms, soln.orbitals, tau1, tau2) -
						soln.atoms[atomnum1].crfg2d(soln.atoms[atomnum2], tau1, tau2);

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
		SimpleMatrix[][] alphabeta;
		if (soln.rm.nonvAlpha + soln.rm.nonvBeta == 0) {
			SimpleMatrix sm = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
			alphabeta = new SimpleMatrix[fockderivstaticalpha.length][2];
			for (int i = 0; i < fockderivstaticalpha.length; i++) {
				alphabeta[i][0] = sm;
				alphabeta[i][1] = sm;
			}
		}
		else alphabeta = Batcher.apply(fockderivstaticalpha, fockderivstaticbeta, SimpleMatrix[][].class,
				State.config.poplethiel_batch_size, (subseta, subsetb) -> {
					SimpleMatrix[] sms = PopleThiel.pt(soln, PopleThiel.toMO(soln.CtaOcc, soln.CaVirt, subseta),
							PopleThiel.toMO(soln.CtbOcc, soln.CbVirt, subsetb));
					SimpleMatrix[][] results = new SimpleMatrix[sms.length][];

					for (int i = 0; i < sms.length; i++) {
						results[i] = PopleThiel.densityDeriv(soln, sms[i]);
					}

					return results;
				});

		int length = alphabeta.length;
		SimpleMatrix hessian = new SimpleMatrix(length, length);

		int[][] indices = new int[alphabeta.length * (alphabeta.length + 1) / 2][];
		int count = 0;
		for (int i = 0; i < alphabeta.length; i++) {
			for (int j = i; j < alphabeta.length; j++) {
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
					E += IntStream.range(0, soln.atoms.length).parallel().mapToDouble(a -> {
						double E2 = 0;

						if (a != atomnum1) {
							E2 += Ederiv2(atomnum1, a, soln.orbsOfAtom, soln.alphaDensity(), soln.betaDensity(),
									soln.atoms, soln.orbitals, tau1, tau2);
							E2 += soln.atoms[atomnum1].crfg2d(soln.atoms[a], tau1, tau2);
						}

						return E2;
					}).sum();
				}
				else {
					E = -Ederiv2(atomnum1, atomnum2, soln.orbsOfAtom, soln.alphaDensity(), soln.betaDensity(),
							soln.atoms, soln.orbitals, tau1, tau2) -
							soln.atoms[atomnum1].crfg2d(soln.atoms[atomnum2], tau1, tau2);
				}

				for (int I = 0; I < soln.orbitals.length; I++) {
					for (int J = 0; J < soln.orbitals.length; J++) {
						E += fockderivstaticalpha[i].get(I, J) * alphabeta[j][0].get(I, J);
						E += fockderivstaticbeta[i].get(I, J) * alphabeta[j][1].get(I, J);
					}
				}

				hessian.set(i, j, E);
				hessian.set(j, i, E);
			}
		});

		return hessian;
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
