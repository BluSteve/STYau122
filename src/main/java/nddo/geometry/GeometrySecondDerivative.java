package nddo.geometry;

import nddo.NDDOAtom;
import nddo.NDDOOrbital;
import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.RecursiveAction;

import static nddo.State.nom;

public class GeometrySecondDerivative {
	private static double Ederiv2(int atomnum1, int atomnum2, int[][] index,
								  SimpleMatrix densityMatrix, NDDOAtom[] atoms,
								  NDDOOrbital[] orbitals, int tau1, int tau2) {

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += densityMatrix.get(i, j) * atoms[atomnum2]
							.Vderiv2(orbitals[i], orbitals[j], tau1, tau2);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e += densityMatrix.get(k, l) * atoms[atomnum1]
							.Vderiv2(orbitals[k], orbitals[l], tau1, tau2);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * densityMatrix.get(i, k) *
							nom.betaderiv2(orbitals[i], orbitals[k], tau1,
									tau2);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {
							e += (densityMatrix.get(i, j) *
									densityMatrix.get(k, l) -
									densityMatrix.get(i, k) * 0.5 *
											densityMatrix.get(j, l))
									* nom
									.getGderiv2(orbitals[i], orbitals[j],
											orbitals[k], orbitals[l], tau1,
											tau2);
						}
					}
				}
			}
		}

		return e;


	}

	private static double Ederiv2(int atomnum1, int atomnum2, int[][] index,
								  SimpleMatrix alphaDensity,
								  SimpleMatrix betaDensity, NDDOAtom[] atoms,
								  NDDOOrbital[] orbitals, int tau1, int tau2) {

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += (alphaDensity.get(i, j) + betaDensity.get(i, j)) *
							atoms[atomnum2]
									.Vderiv2(orbitals[i], orbitals[j], tau1,
											tau2);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e += (alphaDensity.get(k, l) + betaDensity.get(k, l)) *
							atoms[atomnum1]
									.Vderiv2(orbitals[k], orbitals[l], tau1,
											tau2);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * (alphaDensity.get(i, k) + betaDensity.get(i, k)) *
							nom.betaderiv2(orbitals[i], orbitals[k], tau1,
									tau2);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {
							e += ((alphaDensity.get(i, j) +
									betaDensity.get(i, j)) *
									(alphaDensity.get(k, l) +
											betaDensity.get(k, l)) -
									alphaDensity.get(i, k) *
											alphaDensity.get(j, l) -
									betaDensity.get(i, k) *
											betaDensity.get(j, l))
									* nom
									.getGderiv2(orbitals[i], orbitals[j],
											orbitals[k], orbitals[l], tau1,
											tau2);
						}
					}
				}
			}
		}

		return e;


	}

	@Deprecated
	public static double hessianfinite(NDDOAtom[] atoms, SolutionR soln,
									   int atomnum1, int tau1, int atomnum2,
									   int tau2) {

		double initval =
				GeometryDerivative.gradient(atoms, soln, atomnum1, tau1);

		NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum2, tau2);

		SolutionR newsoln = (SolutionR) soln.withNewAtoms(newatoms);

		double finalval =
				GeometryDerivative.gradient(newatoms, newsoln, atomnum1, tau1);

		return 1E7 * (finalval - initval);


	}

	@Deprecated
	public static double hessianfinite(NDDOAtom[] atoms, SolutionU soln,
									   int atomnum1, int tau1, int atomnum2,
									   int tau2) {

		double initval = GeometryDerivative
				.gradientUnrestricted(atoms, soln, atomnum1, tau1);

		NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum2, tau2);

		SolutionU newsoln = (SolutionU) soln.withNewAtoms(newatoms);

		double finalval = GeometryDerivative
				.gradientUnrestricted(newatoms, newsoln, atomnum1, tau1);

		return 1E7 * (finalval - initval);


	}

	public static SimpleMatrix hessianRoutine(SolutionR soln,
											  SimpleMatrix[] fockDerivStatic) {
		SimpleMatrix[] densityDerivs =
				new SimpleMatrix[fockDerivStatic.length];
		int elapsedSize = 0;
		double cores = Runtime.getRuntime().availableProcessors();
		int size = Math.max((int) Math.ceil(fockDerivStatic.length / cores),
				1);

		List<RecursiveAction> subtasks = new ArrayList<>();
		// partitions densityDerivs into batches.
		while (elapsedSize < fockDerivStatic.length) {
			int finalElapsedSize = elapsedSize;
			subtasks.add(new RecursiveAction() {
				@Override
				protected void compute() {
					SimpleMatrix[] subset = Arrays.copyOfRange(fockDerivStatic,
							finalElapsedSize, Math.min(fockDerivStatic.length,
									finalElapsedSize + size));

					SimpleMatrix[] output;
					try {
						output = densityDerivPople(soln, subset);
					} catch (SingularMatrixException e) {
						output = densityDerivThiel(soln, subset);
					}

					// removed .dup() here
					System.arraycopy(output, 0, densityDerivs,
							finalElapsedSize, output.length);
				}
			});
			elapsedSize += size;
		}
		ForkJoinTask.invokeAll(subtasks);

		SimpleMatrix hessian =
				new SimpleMatrix(densityDerivs.length, densityDerivs.length);
		for (int i = 0; i < hessian.numRows(); i++) {
			for (int j = i; j < hessian.numRows(); j++) {
				double E = 0;
				int atomnum1 = i / 3;
				int atomnum2 = j / 3;
				int tau1 = i - 3 * atomnum1;
				int tau2 = j - 3 * atomnum2;

				if (atomnum1 == atomnum2) {
					for (int a = 0; a < soln.atoms.length; a++) {
						if (a != atomnum1) {
							E += Ederiv2(atomnum1, a, soln.orbsOfAtom,
									soln.densityMatrix(), soln.atoms,
									soln.orbitals,
									tau1, tau2);
							E += soln.atoms[atomnum1]
									.crfDeriv2(soln.atoms[a], tau1, tau2);
						}
					}
				}
				else {
					E = -Ederiv2(atomnum1, atomnum2, soln.orbsOfAtom,
							soln.densityMatrix(), soln.atoms, soln.orbitals,
							tau1,
							tau2) - soln.atoms[atomnum1]
							.crfDeriv2(soln.atoms[atomnum2], tau1, tau2);

				}

				for (int I = 0; I < soln.orbitals.length; I++) {
					for (int J = 0; J < soln.orbitals.length; J++) {
						E += fockDerivStatic[i].get(I, J) *
								densityDerivs[j].get(I, J);
					}
				}

				hessian.set(i, j, E);
				hessian.set(j, i, E);
			}
		}
		return hessian;
	}

	public static SimpleMatrix hessianRoutine(SolutionU soln,
											  SimpleMatrix[] fockderivstaticalpha,
											  SimpleMatrix[] fockderivstaticbeta) {


		SimpleMatrix[] densityderivsalpha =
				new SimpleMatrix[fockderivstaticalpha.length];
		SimpleMatrix[] densityderivsbeta =
				new SimpleMatrix[fockderivstaticbeta.length];


		int count = 0;

		for (int a = 0; a < soln.atoms.length; a++) {
			for (int tau = 0; tau < 3; tau++) {

				SimpleMatrix[] matrices = GeometryDerivative
						.densitymatrixderivfinite(soln.atoms, soln, a, tau);
				densityderivsalpha[count] = matrices[0];
				densityderivsbeta[count] = matrices[1];
				count++;
			}
		}

		SimpleMatrix hessian = new SimpleMatrix(densityderivsalpha.length,
				densityderivsalpha.length);

		for (int i = 0; i < hessian.numRows(); i++) {
			for (int j = i; j < hessian.numRows(); j++) {

				double E = 0;

				int atomnum1 = i / 3;

				int atomnum2 = j / 3;

				int tau1 = i - 3 * atomnum1;

				int tau2 = j - 3 * atomnum2;

				if (atomnum1 == atomnum2) {
					for (int a = 0; a < soln.atoms.length; a++) {
						if (a != atomnum1) {
							E += Ederiv2(atomnum1, a, soln.orbsOfAtom,
									soln.alphaDensity(), soln.betaDensity(),
									soln.atoms, soln.orbitals, tau1, tau2);
							E += soln.atoms[atomnum1]
									.crfDeriv2(soln.atoms[a], tau1, tau2);
						}
					}
				}
				else {
					E = -Ederiv2(atomnum1, atomnum2, soln.orbsOfAtom,
							soln.alphaDensity(), soln.betaDensity(),
							soln.atoms,
							soln.orbitals, tau1, tau2) - soln.atoms[atomnum1]
							.crfDeriv2(soln.atoms[atomnum2], tau1, tau2);

				}


				for (int I = 0; I < soln.orbitals.length; I++) {
					for (int J = 0; J < soln.orbitals.length; J++) {
						E += fockderivstaticalpha[i].get(I, J) *
								densityderivsalpha[j].get(I, J);
						E += fockderivstaticbeta[i].get(I, J) *
								densityderivsbeta[j].get(I, J);
					}
				}

				hessian.set(i, j, E);
				hessian.set(j, i, E);
			}
		}

		return hessian;


	}

	public static SimpleMatrix[] densityDerivThiel(SolutionR soln,
												   SimpleMatrix[] fockderivstatic) {
		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix[] xarray = new SimpleMatrix[fockderivstatic.length];
		SimpleMatrix[] rarray = new SimpleMatrix[fockderivstatic.length];
		SimpleMatrix[] dirs = new SimpleMatrix[fockderivstatic.length];

		double[] arrpreconditioner = new double[NOcc * NVirt];

		int counter = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = -soln.E.get(i) - soln.E.get(NOcc + j);

				arrpreconditioner[counter] = Math.pow(e, -0.5);

				counter++;
			}
		}

		SimpleMatrix D = SimpleMatrix.diag(arrpreconditioner);

		for (int a = 0; a < xarray.length; a++) {
			SimpleMatrix F = new SimpleMatrix(NOcc * NVirt, 1);

			int count1 = 0;

			for (int i = 0; i < NOcc; i++) {
				for (int j = 0; j < NVirt; j++) {

					double element = 0;

					for (int u = 0; u < soln.orbitals.length; u++) {
						for (int v = 0; v < soln.orbitals.length; v++) {
							element += soln.C.get(i, u) * soln.C.get(j + NOcc,
									v) * fockderivstatic[a].get(u, v);
						}
					}
					F.set(count1, 0, element);

					count1++;
				}
			}

			F = D.mult(F);

			xarray[a] = new SimpleMatrix(NOcc * NVirt, 1);
			rarray[a] = F;
			dirs[a] = F;
		}


		if (dirs[0].numRows() == 0) {
			SimpleMatrix[] densityderivs =
					new SimpleMatrix[fockderivstatic.length];

			for (int i = 0; i < densityderivs.length; i++) {
				densityderivs[i] = new SimpleMatrix(0, 0);
			}

			return densityderivs;
		}


		while (Utils.numNotNull(rarray) > 0) {
			ArrayList<SimpleMatrix> d = new ArrayList<>();

			ArrayList<SimpleMatrix> p = new ArrayList<>();

			for (int i = 0; i < rarray.length; i++) {
				if (rarray[i] != null) {
					d.add(new SimpleMatrix(dirs[i]));
					p.add(D.mult(
							computeResponseVectorsThiel(dirs[i], soln)));
				}
			}

			SimpleMatrix solver =
					new SimpleMatrix(p.size(), p.size());
			SimpleMatrix rhsvec =
					new SimpleMatrix(p.size(), rarray.length);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					double[] arrrhs = new double[p.size()];

					for (int i = 0; i < arrrhs.length; i++) {
						arrrhs[i] = 2 *
								rarray[a].transpose().mult(d.get(i))
										.get(0, 0);

					}
					rhsvec.setColumn(a, 0, arrrhs);
				}
			}

			for (int i = 0; i < solver.numRows(); i++) {
				for (int j = i; j < solver.numRows(); j++) {
					double val2 =
							p.get(j).transpose().mult(d.get(i))
									.get(0, 0) + p.get(i).transpose()
									.mult(d.get(j)).get(0, 0);

					solver.set(i, j, val2);
					solver.set(j, i, val2);
				}
			}

			SimpleMatrix alpha;
			try {
				alpha = solver.solve(rhsvec);
			} catch (SingularMatrixException e) {
				alpha = SimpleMatrix.ones(solver.numCols(), rhsvec.numCols());
			}

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					for (int i = 0; i < alpha.numRows(); i++) {
						xarray[a] = xarray[a].plus(
								d.get(i).scale(alpha.get(i, a)));

						rarray[a] = rarray[a].minus(
								p.get(i).scale(alpha.get(i, a)));

					}

					if (Utils.mag(rarray[a]) < 1E-6) {//todo change this if you want
						rarray[a] = null;
					}
					else {
//						System.out.println("convergence test: " + mag
//								(rarray[a]));
					}
				}
			}

			solver = new SimpleMatrix(solver.numRows(),
					solver.numRows());

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					double[] arrrhs = new double[solver.numRows()];

					for (int i = 0; i < arrrhs.length; i++) {
						arrrhs[i] = -rarray[a].transpose()
								.mult(p.get(i)).get(0, 0);

					}
					rhsvec.setColumn(a, 0, arrrhs);
				}
			}

			for (int i = 0; i < solver.numRows(); i++) {
				for (int j = 0; j < solver.numRows(); j++) {
					solver.set(i, j,
							d.get(j).transpose().mult(p.get(i))
									.get(0, 0));
				}
			}

			SimpleMatrix beta = solver.solve(rhsvec);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					dirs[a] = rarray[a];

					for (int i = 0; i < beta.numRows(); i++) {
						dirs[a] = dirs[a].plus(
								d.get(i).scale(beta.get(i, a)));
					}
				}
			}
		}

		SimpleMatrix[] densityMatrixDerivs =
				new SimpleMatrix[fockderivstatic.length];

		for (int a = 0; a < fockderivstatic.length; a++) {
			SimpleMatrix densityMatrixDeriv = new SimpleMatrix
					(soln.orbitals.length, soln.orbitals.length);

			for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
				for (int v = u; v < densityMatrixDeriv.numCols(); v++) {
					double sum = 0;
					int count = 0;
					for (int i = 0; i < NOcc; i++) {
						for (int j = 0; j < NVirt; j++) {
							sum -= 2 * (soln.C.get(i, u) *
									soln.C.get(j + NOcc, v) +
									soln.C.get(j + NOcc, u) *
											soln.C.get(i, v)) *
									xarray[a].get(count, 0);
							count++;
						}
					}

					densityMatrixDeriv.set(u, v, sum);
					densityMatrixDeriv.set(v, u, sum);
				}
			}

			densityMatrixDerivs[a] = densityMatrixDeriv;
		}

		return densityMatrixDerivs;
	}

	public static SimpleMatrix[] getxarrayPople(SolutionR soln,
												SimpleMatrix[] fockderivstatic) {
		int NOcc = (int) (soln.nElectrons / 2.0);
		int NVirt = soln.orbitals.length - NOcc;
		int length = fockderivstatic.length;
		int nonv = NOcc * NVirt;

		if (nonv == 0) {
			SimpleMatrix[] xarray = new SimpleMatrix[length];

			for (int i = 0; i < xarray.length; i++) {
				xarray[i] = new SimpleMatrix(0, 0);
			}

			return xarray;
		}

		// array initialization
		SimpleMatrix[] xarray = new SimpleMatrix[length];
		SimpleMatrix[] barray = new SimpleMatrix[length];
		SimpleMatrix[] parray = new SimpleMatrix[length];
		SimpleMatrix[] Farray = new SimpleMatrix[length];
		SimpleMatrix[] rarray = new SimpleMatrix[length];

		// configure preconditioners
		double[] Darr = new double[nonv];
		double[] Dinvarr = new double[nonv];

		int counter = 0;
		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = -soln.E.get(i) + soln.E.get(NOcc + j);

				Darr[counter] = Math.pow(e, -0.5);
				Dinvarr[counter] = Math.pow(e, 0.5);

				counter++;
			}
		}

		// convert AO to MO basis
		SimpleMatrix F = new SimpleMatrix(nonv, length);
		for (int a = 0; a < length; a++) {
			SimpleMatrix f = new SimpleMatrix(nonv, 1);

			int count = 0;

			for (int i = 0; i < NOcc; i++) { // kappa
				for (int j = 0; j < NVirt; j++) { // i
					double element = 0;

					for (int u = 0; u < soln.orbitals.length; u++) {
						for (int v = 0; v < soln.orbitals.length; v++) {
							element += soln.C.get(i, u) *
									soln.C.get(j + NOcc, v) *
									fockderivstatic[a].get(u, v);
						}
					}

					element /= soln.E.get(j + NOcc) - soln.E.get(i);

					f.set(count, 0, element);

					count++;
				}
			}

			CommonOps_DDRM.multRows(Darr, f.getDDRM());
			barray[a] = f;
			Farray[a] = barray[a].copy();
			F.setColumn(a, 0, barray[a].getDDRM().data);
		}

		// main loop
		int[] iterable = new int[length];

		// 0: B, 1: Bt, 2: Bn, 3: P, 4: BmP
		List<SimpleMatrix[]> prevs = new ArrayList<>();
		List<Double> dots = new ArrayList<>();

		while (Utils.numIterable(iterable) > 0) {
			// orthogonalize barray
			for (int i = 0; i < barray.length; i++) {
				for (int j = 0; j < i; j++) {
					barray[i].plusi(barray[i].dot(barray[j]) /
							barray[j].dot(barray[j]), barray[j].negativei());
				}
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix[] prev = new SimpleMatrix[5];
				prev[0] = barray[i]; // original barray object here
				prev[1] = barray[i].transpose();
				prev[2] = barray[i].negative();
				dots.add(barray[i].dot(barray[i]));

				// parray[i] stays the same object throughout
				SimpleMatrix bc = barray[i].copy();
				CommonOps_DDRM.multRows(Dinvarr, bc.getDDRM());
				SimpleMatrix crv = computeResponseVectorsPople(bc, soln);
				CommonOps_DDRM.multRows(Darr, crv.getDDRM());
				parray[i] = crv;

				prev[3] = parray[i];
				prev[4] = barray[i].minus(parray[i]);

				prevs.add(prev);
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix newb = parray[i].copy();

				// orthogonalize against all previous Bs
				for (int j = 0; j < prevs.size(); j++) {
					SimpleMatrix[] prev = prevs.get(j);
					SimpleMatrix transpose = prev[1];
					double num = transpose.mult(parray[i]).get(0) /
							dots.get(j);

					newb.plusi(num, prev[2]);
				}

				barray[i] = newb; // new barray object created
			}

			SimpleMatrix Bt = new SimpleMatrix(prevs.size(), nonv);
			SimpleMatrix P = new SimpleMatrix(nonv, prevs.size());

			for (int i = 0; i < prevs.size(); i++) {
				Bt.setRow(i, 0, prevs.get(i)[0].getDDRM().data);
				P.setColumn(i, 0, prevs.get(i)[4].getDDRM().data);
			}

			SimpleMatrix rhs = Bt.mult(F);
			SimpleMatrix lhs = Bt.mult(P);
			// alpha dimensions are prevBs x length
			SimpleMatrix alpha = lhs.solve(rhs);

			// reset r and x array
			for (int a = 0; a < length; a++) {
				rarray[a] = new SimpleMatrix(nonv, 1);
				xarray[a] = new SimpleMatrix(nonv, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {
					// B with tilde
					rarray[j].plusi(alpha.get(i, j), prevs.get(i)[4]);
					xarray[j].plusi(alpha.get(i, j), prevs.get(i)[0]);
				}
			}

			for (int j = 0; j < alpha.numCols(); j++) {
				// B0 is Farray, no tilde
				rarray[j].minusi(Farray[j]);
				CommonOps_DDRM.multRows(Dinvarr, xarray[j].getDDRM());

				double rMag = Utils.mag(rarray[j]);
				if (rMag < 1E-7) {
					iterable[j] = 1;
				}
				else if (Double.isNaN(rMag)) {
					soln.getRm().getLogger()
							.warn("Pople algorithm fails; reverting to " +
									"Thiel algorithm (don't panic)...");
					throw new SingularMatrixException();
				}
				else {
					iterable[j] = 0;
				}
			}
		}

		return xarray;
	}

	public static SimpleMatrix[] densityDerivPople(SolutionR soln,
												   SimpleMatrix[] fockderivstatic) {
		int NOcc = (int) (soln.nElectrons / 2.0);
		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix[] xarray = getxarrayPople(soln, fockderivstatic);

		SimpleMatrix[] densityMatrixDerivs =
				new SimpleMatrix[fockderivstatic.length];

		for (int a = 0; a < fockderivstatic.length; a++) {
			SimpleMatrix densityMatrixDeriv =
					new SimpleMatrix(soln.orbitals.length,
							soln.orbitals.length);

			for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
				for (int v = u; v < densityMatrixDeriv.numCols(); v++) {
					double sum = 0;
					int count = 0;
					for (int i = 0; i < NOcc; i++) {
						for (int j = 0; j < NVirt; j++) {
							sum -= 2 * (soln.C.get(i, u) *
									soln.C.get(j + NOcc, v) +
									soln.C.get(j + NOcc, u) *
											soln.C.get(i, v)) *
									xarray[a].get(count, 0);
							count++;
						}
					}

					densityMatrixDeriv.set(u, v, sum);
					densityMatrixDeriv.set(v, u, sum);
				}
			}

			densityMatrixDerivs[a] = densityMatrixDeriv;
		}

		return densityMatrixDerivs;
	}


	private static SimpleMatrix computeResponseVectorsThiel(SimpleMatrix x,
															SolutionR soln) {

		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix densityMatrixDeriv =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
			for (int v = u; v < densityMatrixDeriv.numCols(); v++) {
				double sum = 0;
				int count = 0;
				for (int i = 0; i < NOcc; i++) {
					for (int j = 0; j < NVirt; j++) {
						sum -= 2 * (soln.C.get(i, u) * soln.C.get(j + NOcc,
								v) +
								soln.C.get(j + NOcc, u) * soln.C.get(i, v)) *
								x.get(count, 0);
						count++;
					}
				}

				densityMatrixDeriv.set(u, v, sum);
				densityMatrixDeriv.set(v, u, sum);
			}
		}


		SimpleMatrix responsematrix =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = j; k < soln.orbitals.length; k++) {
				double val = 0;
				if (j == k) {

					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							val += densityMatrixDeriv.get(l, l) *
									integralArray[integralcount];
							integralcount++;
						}
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] ==
											soln.atomOfOrb[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += densityMatrixDeriv.get(j, k) *
							integralArray[integralcount];
					integralcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] ==
											soln.atomOfOrb[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.orbsOfAtom[soln.atomOfOrb[k]]) {
								if (m > -1) {
									val += densityMatrixDeriv.get(l, m) *
											integralArray[integralcount];
									integralcount++;
								}
							}
						}
					}
				}

				responsematrix.set(j, k, val);
				responsematrix.set(k, j, val);
			}
		}

		SimpleMatrix R = new SimpleMatrix(NOcc * NVirt, 1);

		int count1 = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {

				double element = 0;

				for (int u = 0; u < soln.orbitals.length; u++) {
					for (int v = 0; v < soln.orbitals.length; v++) {
						element += soln.C.get(i, u) * soln.C.get(j + NOcc, v) *
								responsematrix.get(u, v);
					}
				}


				R.set(count1, 0, element);

				count1++;
			}
		}


		SimpleMatrix p = new SimpleMatrix(NOcc * NVirt, 1);

		int counter = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				p.set(counter, 0, -R.get(counter, 0) +
						(soln.E.get(j + NOcc) - soln.E.get(i)) *
								x.get(counter));
				counter++;
			}
		}


		return p;
	}

	public static SimpleMatrix computeResponseVectorsPople(SimpleMatrix x,
														   SolutionR soln) {
		int NOcc = (int) (soln.nElectrons / 2.0);

		int NVirt = soln.orbitals.length - NOcc;

		SimpleMatrix densityMatrixDeriv =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		for (int u = 0; u < densityMatrixDeriv.numRows(); u++) {
			for (int v = u; v < densityMatrixDeriv.numCols(); v++) {
				double sum = 0;
				int count = 0;
				for (int i = 0; i < NOcc; i++) {
					for (int j = 0; j < NVirt; j++) {
						sum -= 2 * (soln.C.get(i, u) * soln.C.get(j + NOcc,
								v) +
								soln.C.get(j + NOcc, u) * soln.C.get(i, v)) *
								x.get(count, 0);
						count++;
					}
				}

				densityMatrixDeriv.set(u, v, sum);
				densityMatrixDeriv.set(v, u, sum);
			}
		}


		SimpleMatrix responsematrix =
				new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = j; k < soln.orbitals.length; k++) {
				double val = 0;
				if (j == k) {

					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							val += densityMatrixDeriv.get(l, l) *
									integralArray[integralcount];
							integralcount++;
						}
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] ==
											soln.atomOfOrb[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += densityMatrixDeriv.get(j, k) *
							integralArray[integralcount];
					integralcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.missingOfAtom[soln.atomOfOrb[j]]) {
								if (m > -1) {
									if (soln.atomOfOrb[l] ==
											soln.atomOfOrb[m]) {
										val += densityMatrixDeriv.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}

							}
						}
					}
				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						if (l > -1) {
							for (int m :
									soln.orbsOfAtom[soln.atomOfOrb[k]]) {
								if (m > -1) {
									val += densityMatrixDeriv.get(l, m) *
											integralArray[integralcount];
									integralcount++;
								}
							}
						}
					}
				}

				responsematrix.set(j, k, val);
				responsematrix.set(k, j, val);
			}
		}

		SimpleMatrix R = new SimpleMatrix(NOcc * NVirt, 1);

		int count1 = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {

				double element = 0;

				for (int u = 0; u < soln.orbitals.length; u++) {
					for (int v = 0; v < soln.orbitals.length; v++) {
						element += soln.C.get(i, u) * soln.C.get(j + NOcc, v) *
								responsematrix.get(u, v);
					}
				}

				element = element / (soln.E.get(j + NOcc) - soln.E.get(i));


				R.set(count1, 0, element);

				count1++;
			}
		}


		return R;
	}
}
