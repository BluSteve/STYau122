package nddo.geometry;

import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.math.PopleThiel;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.apache.commons.lang3.time.StopWatch;
import org.ejml.simple.SimpleMatrix;
import tools.Batcher;
import tools.Pow;
import tools.Utils;

import java.util.ArrayList;

import static nddo.State.nom;
import static tools.Utils.mag;

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
		else densityDerivs = Batcher.apply(fockderivstatic, subset -> densityDeriv(soln, subset));

		SimpleMatrix hessian = new SimpleMatrix(densityDerivs.length, densityDerivs.length);

		int[][] indices = new int[hessian.numRows() * (hessian.numRows() + 1) / 2][];
		int count = 0;
		for (int i = 0; i < hessian.numRows(); i++) {
			for (int j = i; j < hessian.numCols(); j++) {
				indices[count] = new int[]{i, j};
				count++;
			}
		}

		Batcher.consume(indices, subset -> {
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
		else sms = Batcher.apply(new SimpleMatrix[][]{fockderivstaticalpha, fockderivstaticbeta},
				subset -> densityDeriv(soln, subset));

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

		Batcher.consume(indices, subset -> {
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

	private static SimpleMatrix[][] densityDeriv(SolutionU soln, SimpleMatrix[][] fockderivstatics) {
		SimpleMatrix[] fockderivstaticalpha = fockderivstatics[0];
		SimpleMatrix[] fockderivstaticbeta = fockderivstatics[1];

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


	@Deprecated
	private static SimpleMatrix[][] densityDerivPople(SolutionU soln, SimpleMatrix[] fockderivstaticalpha,
													  SimpleMatrix[] fockderivstaticbeta) {

		StopWatch sw = new StopWatch();
		sw.start();
		int NOccAlpha = soln.rm.nOccAlpha;
		int NOccBeta = soln.rm.nOccBeta;

		int NVirtAlpha = soln.rm.nVirtAlpha;
		int NVirtBeta = soln.rm.nVirtBeta;

		SimpleMatrix[] xarray = new SimpleMatrix[fockderivstaticalpha.length];
		SimpleMatrix[] barray = new SimpleMatrix[fockderivstaticalpha.length];
		SimpleMatrix[] parray = new SimpleMatrix[fockderivstaticalpha.length];
		SimpleMatrix[] Farray = new SimpleMatrix[fockderivstaticalpha.length];
		SimpleMatrix[] rarray = new SimpleMatrix[fockderivstaticalpha.length];

		SimpleMatrix preconditioner = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);
		SimpleMatrix preconditionerinv = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);

		int counter = 0;

		for (int i = 0; i < NOccAlpha; i++) {
			for (int j = 0; j < NVirtAlpha; j++) {
				double e = (-soln.Ea.get(i) + soln.Ea.get(NOccAlpha + j));
				preconditioner.set(counter, Pow.pow(e, -0.5));
				preconditionerinv.set(counter, Pow.pow(e, 0.5));
				counter++;
			}
		}

		for (int i = 0; i < NOccBeta; i++) {
			for (int j = 0; j < NVirtBeta; j++) {
				double e = (-soln.Eb.get(i) + soln.Eb.get(NOccBeta + j));
				preconditioner.set(counter, Pow.pow(e, -0.5));
				preconditionerinv.set(counter, Pow.pow(e, 0.5));
				counter++;
			}
		}

		final SimpleMatrix D = SimpleMatrix.diag(preconditioner.getDDRM().data);

		final SimpleMatrix Dinv = SimpleMatrix.diag(preconditionerinv.getDDRM().data);

//      SimpleMatrix D = SimpleMatrix.eye(NOcc * NVirt);
//
//      SimpleMatrix Dinv = SimpleMatrix.eye(NOcc * NVirt);

		for (int a = 0; a < xarray.length; a++) {
			SimpleMatrix F = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);

			int count1 = 0;

			for (int i = 0; i < NOccAlpha; i++) { // kappa
				for (int j = 0; j < NVirtAlpha; j++) { // i

					double element = 0;

					for (int u = 0; u < soln.orbitals.length; u++) {
						for (int v = 0; v < soln.orbitals.length; v++) {
							element += soln.Cta.get(i, u) * soln.Cta.get(j + NOccAlpha, v) *
									fockderivstaticalpha[a].get(u, v);
						}
					}

					element = element / (soln.Ea.get(j + NOccAlpha) - soln.Ea.get(i));


					F.set(count1, 0, element);

					count1++;
				}
			}

			for (int i = 0; i < NOccBeta; i++) { // kappa
				for (int j = 0; j < NVirtBeta; j++) { // i

					double element = 0;

					for (int u = 0; u < soln.orbitals.length; u++) {
						for (int v = 0; v < soln.orbitals.length; v++) {
							element +=
									soln.Ctb.get(i, u) * soln.Ctb.get(j + NOccBeta, v) * fockderivstaticbeta[a].get(u,
											v);
						}
					}

					element = element / (soln.Eb.get(j + NOccBeta) - soln.Eb.get(i));


					F.set(count1, 0, element);

					count1++;
				}
			}

			F = D.mult(F);

			xarray[a] = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);
			rarray[a] = xarray[a].copy();
			barray[a] = F.copy();
			Farray[a] = F.copy();
		}


		if (barray[0].numRows() == 0) {
			SimpleMatrix[] densityderivs = new SimpleMatrix[fockderivstaticalpha.length];

			for (int i = 0; i < densityderivs.length; i++) {
				densityderivs[i] = new SimpleMatrix(0, 0);
			}

			return new SimpleMatrix[][]{densityderivs, densityderivs};
		}


		ArrayList<SimpleMatrix> prevBs = new ArrayList<>();

		ArrayList<SimpleMatrix> prevPs = new ArrayList<>();

		int[] iterable = new int[barray.length];


		SimpleMatrix F = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, Farray.length);


		for (int i = 0; i < Farray.length; i++) {
			F.setColumn(i, 0, Farray[i].getDDRM().data);
		}

		int numit = 0;

		while (Utils.numIterable(iterable) > 0) {

			numit++;

			for (int i = 0; i < barray.length; i++) {
				for (int j = 0; j < i; j++) {
					barray[i] = barray[i].minus(barray[j]
							.scale(barray[i].dot(barray[j]) /
									barray[j].dot(barray[j])));
				}
			}

			System.out.println("Geom only " + Utils.numIterable(iterable) + " left to go!");

			for (int i = 0; i < barray.length; i++) {

				prevBs.add(barray[i].copy());
				parray[i] = D.mult(computeResponseVectorsPople(Dinv.mult(barray[i].copy()), soln));
				prevPs.add(parray[i].copy());
			}

			for (int i = 0; i < barray.length; i++) {

				SimpleMatrix newb = parray[i];

				for (SimpleMatrix prevB : prevBs) {
					double num = prevB.transpose().mult(parray[i]).get(0) / prevB.transpose().mult(prevB).get(0);

					newb = newb.minus(prevB.scale(num));
				}

				barray[i] = newb.copy();


			}

			SimpleMatrix B = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, prevBs.size());
			SimpleMatrix P = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, prevBs.size());

			for (int i = 0; i < prevBs.size(); i++) {

				B.setColumn(i, 0, prevBs.get(i).getDDRM().data);

				P.setColumn(i, 0, prevBs.get(i).minus(prevPs.get(i)).getDDRM().data);

			}


			SimpleMatrix lhs = B.transpose().mult(P);

			SimpleMatrix rhs = B.transpose().mult(F);

			SimpleMatrix alpha = lhs.solve(rhs);

			for (int a = 0; a < xarray.length; a++) {

				rarray[a] = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);
				xarray[a] = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {
					// B with tilde
					rarray[j] =
							rarray[j].plus((prevBs.get(i).minus(prevPs.get(i)))
									.scale(alpha.get(i, j)));
					xarray[j] =
							xarray[j].plus(prevBs.get(i).scale(alpha.get(i, j)));
				}
			}

			for (int j = 0; j < alpha.numCols(); j++) {

				// B0 is Farray, no tilde
				rarray[j] = rarray[j].minus(Farray[j]);

				xarray[j] = Dinv.mult(xarray[j]);

				if (mag(rarray[j]) < 1E-8) {
					iterable[j] = 1;
				}
				else if (Double.isNaN(mag(rarray[j]))) {
					System.err.println(
							"Pople algorithm fails; reverting to Thiel " +
									"algorithm (don't panic)...");
					return densityDeriv(soln, new SimpleMatrix[][]{fockderivstaticalpha, fockderivstaticbeta});
				}
				else {
					iterable[j] = 0;

				}
				System.out.println("convergence test: " + mag(rarray[j]));

			}


		}


		SimpleMatrix[] densityderivsalpha =
				new SimpleMatrix[fockderivstaticalpha.length];

		SimpleMatrix[] densityderivsbeta =
				new SimpleMatrix[fockderivstaticbeta.length];


		for (int a = 0; a < fockderivstaticalpha.length; a++) {

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
									xarray[a].get(count, 0);
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
									xarray[a].get(count, 0);
							count++;
						}
					}

					densityderivbeta.set(u, v, sum);
					densityderivbeta.set(v, u, sum);
				}
			}

			densityderivsalpha[a] = densityderivalpha;
			densityderivsbeta[a] = densityderivbeta;

		}

		System.out.println("numit (Pople): " + numit);

		//System.exit(0);

//		System.err.println("Time: " + sw.getTime());


		return new SimpleMatrix[][]{densityderivsalpha, densityderivsbeta};


	}

	@Deprecated
	private static SimpleMatrix computeResponseVectorsPople(SimpleMatrix xarray, SolutionU soln) {

		int NOccAlpha = soln.rm.nOccAlpha;
		int NOccBeta = soln.rm.nOccBeta;

		int NVirtAlpha = soln.rm.nVirtAlpha;
		int NVirtBeta = soln.rm.nVirtBeta;

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
								xarray.get(count, 0);
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
								xarray.get(count, 0);
						count++;
					}
				}

				densityderivbeta.set(u, v, sum);
				densityderivbeta.set(v, u, sum);
			}
		}

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


		SimpleMatrix R = new SimpleMatrix(NOccAlpha * NVirtAlpha + NOccBeta * NVirtBeta, 1);

		int count1 = 0;

		for (int i = 0; i < NOccAlpha; i++) {
			for (int j = 0; j < NVirtAlpha; j++) {

				double element = 0;

				for (int u = 0; u < soln.orbitals.length; u++) {
					for (int v = 0; v < soln.orbitals.length; v++) {
						element += soln.Cta.get(i, u) * soln.Cta.get(j + NOccAlpha, v) *
								responsealpha.get(u, v);
					}
				}

				element = element / (soln.Ea.get(j + NOccAlpha) - soln.Ea.get(i));


				R.set(count1, 0, element);

				count1++;
			}
		}

		for (int i = 0; i < NOccBeta; i++) {
			for (int j = 0; j < NVirtBeta; j++) {

				double element = 0;

				for (int u = 0; u < soln.orbitals.length; u++) {
					for (int v = 0; v < soln.orbitals.length; v++) {
						element += soln.Ctb.get(i, u) * soln.Ctb.get(j + NOccBeta, v) *
								responsebeta.get(u, v);
					}
				}

				element = element / (soln.Eb.get(j + NOccBeta) - soln.Eb.get(i));


				R.set(count1, 0, element);

				count1++;
			}
		}


		return R;
	}
}
