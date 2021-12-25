package nddo.math;

import nddo.solution.SolutionR;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.SpecializedOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.ejml.dense.row.CommonOps_DDRM.multRows;
import static tools.Utils.mag;

public class PopleThiel {
	public static SimpleMatrix[] pople(SolutionR soln, SimpleMatrix[] fockderivstatic) {
		// Pople alg will solve any equation of the form (1-D)x=B
		int NOcc = soln.rm.nOccAlpha;
		int NVirt = soln.rm.nVirtAlpha;
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
		SimpleMatrix[] barray = new SimpleMatrix[length]; // B tilde
		SimpleMatrix[] parray = new SimpleMatrix[length]; // trial R / (ej - ei) as well as D * B tilde
		SimpleMatrix[] Farray = new SimpleMatrix[length]; // F in MO basis divided by (ej - ei)
		SimpleMatrix[] rarray = new SimpleMatrix[length]; // x - F
		double[] oldrMags = new double[rarray.length];
		Arrays.fill(oldrMags, 1);

		// configure preconditioners
		double[] Darr = new double[nonv];
		double[] Dinvarr = new double[nonv];

		int counter = 0;
		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = -soln.E.get(i) + soln.E.get(NOcc + j);

				Dinvarr[counter] = Math.sqrt(e);
				Darr[counter] = 1 / Dinvarr[counter];

				counter++;
			}
		}

		SimpleMatrix F = new SimpleMatrix(nonv, length); // scaled fockderivstatic vectors in matrix form, cf. Farray
		for (int a = 0; a < length; a++) {
			SimpleMatrix Foccvirt = soln.CtOcc.mult(fockderivstatic[a]).mult(soln.CVirt); // convert AO to MO basis
			SimpleMatrix f = Foccvirt.elementDivi(soln.Emat); // divided by (ej - ei)
			f.reshape(nonv, 1);

			multRows(Darr, f.getDDRM());
			barray[a] = f;
			Farray[a] = barray[a].copy();
			F.setColumn(a, 0, barray[a].getDDRM().data);
		}

		// main loop
		boolean[] iterable = new boolean[length];
		boolean[] looselyIterable = new boolean[length];

		// 0: B, 1: Bt, 2: -B, 3: P, 4: B-P
		List<SimpleMatrix[]> prevs = new ArrayList<>();
		List<Double> dots = new ArrayList<>();

		SimpleMatrix alpha = null;

		bigLoop:
		while (numIterable(iterable) > 0) {
			// orthogonalize barray
			for (int i = 1; i < barray.length; i++) {
				for (int j = 0; j < i; j++) {
					barray[i].plusi(barray[i].dot(barray[j]) / barray[j].dot(barray[j]), barray[j].negative());
				}
			}

			for (int i = 0; i < length; i++) {
				// parray[i] stays the same object throughout
				SimpleMatrix b = barray[i].copy();
				multRows(Dinvarr, b.getDDRM());
				SimpleMatrix p = computeResponseVectorsPople(soln, b); // P = D * B
				multRows(Darr, p.getDDRM());
				parray[i] = p;

				SimpleMatrix[] prev = new SimpleMatrix[5];
				prev[0] = barray[i]; // original barray object here
				prev[1] = barray[i].transpose();
				prev[2] = barray[i].negative();
				dots.add(barray[i].dot(barray[i]));

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
					double num = transpose.mult(parray[i]).get(0) / dots.get(j);

					newb.plusi(num, prev[2]);
				}

				barray[i] = newb; // new barray object created
			}

			SimpleMatrix Bt = new SimpleMatrix(prevs.size(), nonv);
			SimpleMatrix BminusP = new SimpleMatrix(nonv, prevs.size());

			for (int i = 0; i < prevs.size(); i++) {
				Bt.setRow(i, 0, prevs.get(i)[0].getDDRM().data);
				BminusP.setColumn(i, 0, prevs.get(i)[4].getDDRM().data);
			}

			SimpleMatrix rhs = Bt.mult(F);
			SimpleMatrix lhs = Bt.mult(BminusP);
			// alpha dimensions are prevBs x length
			alpha = lhs.solve(rhs);

			// reset r array
			for (int a = 0; a < length; a++) {
				rarray[a] = new SimpleMatrix(nonv, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {
					rarray[j].plusi(alpha.get(i, j), prevs.get(i)[4]);
				}
			}

			for (int j = 0; j < alpha.numCols(); j++) {
				double mag = SpecializedOps_DDRM.diffNormF_fast(rarray[j].getDDRM(), Farray[j].getDDRM());

				if (mag > oldrMags[j] || mag != mag) { // unstable
					if (numIterable(looselyIterable) == 0) {
						soln.getRm().getLogger().warn("Slight numerical instability detected; " +
								"returning lower precision values. rarray mag = {}", mag);
						break bigLoop; // i.e. mag is tolerable
					}

					soln.getRm().getLogger().warn("Pople algorithm fails; reverting to Thiel algorithm...");
					throw new SingularMatrixException();
				}

				if (mag < 1E-7) {
					if (mag < 1E-10) looselyIterable[j] = true;
					oldrMags[j] = mag;
					iterable[j] = true;
				}
				else {
					iterable[j] = false;
				}
			}
		}

		for (int i = 0; i < length; i++) {
			xarray[i] = new SimpleMatrix(nonv, 1);
		}

		for (int i = 0; i < alpha.numRows(); i++) {
			for (int j = 0; j < alpha.numCols(); j++) {
				xarray[j].plusi(alpha.get(i, j), prevs.get(i)[0]);
			}
		}

		for (int j = 0; j < alpha.numCols(); j++) {
			multRows(Dinvarr, xarray[j].getDDRM());
		}

		return xarray;
	}

	public static SimpleMatrix[] thiel(SolutionR soln, SimpleMatrix[] fockderivstatic) {
		int NOcc = soln.rm.nOccAlpha;
		int NVirt = soln.rm.nVirtAlpha;
		int length = fockderivstatic.length;
		int nonv = NOcc * NVirt;

		if (nonv == 0) {
			SimpleMatrix[] densityderivs = new SimpleMatrix[length];

			for (int i = 0; i < densityderivs.length; i++) {
				densityderivs[i] = new SimpleMatrix(0, 0);
			}

			return densityderivs;
		}

		SimpleMatrix[] xarray = new SimpleMatrix[length];
		SimpleMatrix[] rarray = new SimpleMatrix[length];
		SimpleMatrix[] dirs = new SimpleMatrix[length];

		double[] Darr = new double[nonv];

		int counter = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = -soln.E.get(i) - soln.E.get(NOcc + j);

				Darr[counter] = 1 / Math.sqrt(e);

				counter++;
			}
		}

		for (int a = 0; a < xarray.length; a++) {
			SimpleMatrix f = soln.CtOcc.mult(fockderivstatic[a]).mult(soln.CVirt); // convert AO to MO basis
			f.reshape(nonv, 1);
			multRows(Darr, f.getDDRM());

			xarray[a] = new SimpleMatrix(NOcc * NVirt, 1);
			rarray[a] = f;
			dirs[a] = f;
		}

		while (Utils.numNotNull(rarray) > 0) {
			ArrayList<SimpleMatrix> d = new ArrayList<>();
			ArrayList<SimpleMatrix> p = new ArrayList<>();

			for (int i = 0; i < rarray.length; i++) {
				if (rarray[i] != null) {
					d.add(new SimpleMatrix(dirs[i]));

					SimpleMatrix sm = computeResponseVectorsThiel(soln, dirs[i]);
					multRows(Darr, sm.getDDRM());
					p.add(sm);
				}
			}

			SimpleMatrix solver = new SimpleMatrix(p.size(), p.size());
			SimpleMatrix rhsvec = new SimpleMatrix(p.size(), rarray.length);

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					double[] arrrhs = new double[p.size()];

					for (int i = 0; i < arrrhs.length; i++) {
						arrrhs[i] = 2 * rarray[a].transpose().mult(d.get(i)).get(0, 0);

					}
					rhsvec.setColumn(a, 0, arrrhs);
				}
			}

			for (int i = 0; i < solver.numRows(); i++) {
				for (int j = i; j < solver.numRows(); j++) {
					double val2 = p.get(j).transpose().mult(d.get(i)).get(0, 0) +
							p.get(i).transpose().mult(d.get(j)).get(0, 0);

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
						xarray[a] = xarray[a].plus(d.get(i).scale(alpha.get(i, a)));
						rarray[a] = rarray[a].minus(p.get(i).scale(alpha.get(i, a)));
					}

					if (mag(rarray[a]) < 1E-6) { // todo change this if you want
						rarray[a] = null;
					}
				}
			}

			solver = new SimpleMatrix(solver.numRows(), solver.numRows());

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					double[] arrrhs = new double[solver.numRows()];

					for (int i = 0; i < arrrhs.length; i++) {
						arrrhs[i] = -rarray[a].transpose().mult(p.get(i)).get(0, 0);

					}
					rhsvec.setColumn(a, 0, arrrhs);
				}
			}

			for (int i = 0; i < solver.numRows(); i++) {
				for (int j = 0; j < solver.numRows(); j++) {
					solver.set(i, j, d.get(j).transpose().mult(p.get(i)).get(0, 0));
				}
			}

			SimpleMatrix beta;
			try {
				beta = solver.solve(rhsvec);
			} catch (SingularMatrixException e) {
				beta = SimpleMatrix.ones(solver.numCols(), rhsvec.numCols());
			}

			for (int a = 0; a < rhsvec.numCols(); a++) {
				if (rarray[a] != null) {
					dirs[a] = rarray[a];

					for (int i = 0; i < beta.numRows(); i++) {
						dirs[a] = dirs[a].plus(d.get(i).scale(beta.get(i, a)));
					}
				}
			}
		}

		return xarray;
	}

	public static SimpleMatrix computeResponseVectorsPople(SolutionR soln, SimpleMatrix x) {
		// x is B tilde, i.e. a guess
		// p = d * b. this evaluates d * b without finding out what d is, as the latter is slow.
		int NOcc = soln.rm.nOccAlpha;
		int NVirt = soln.rm.nVirtAlpha;

		SimpleMatrix xmat = x.copy();
		xmat.reshape(NOcc, NVirt);

		SimpleMatrix densityMatrixDeriv = soln.COcc.mult(xmat).mult(soln.CtVirt)
				.plusi(soln.CVirt.mult(xmat.transpose().mult(soln.CtOcc))).scalei(-2);

		SimpleMatrix responsematrix = new SimpleMatrix(soln.orbitals.length, soln.orbitals.length);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.orbitals.length; j++) {
			for (int k = j; k < soln.orbitals.length; k++) {
				double val = 0;
				if (j == k) {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						val += densityMatrixDeriv.get(l, l) * integralArray[integralcount];
						integralcount++;
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += densityMatrixDeriv.get(l, m) * integralArray[integralcount];
								integralcount++;
							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += densityMatrixDeriv.get(j, k) * integralArray[integralcount];
					integralcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += densityMatrixDeriv.get(l, m) * integralArray[integralcount];
								integralcount++;
							}
						}
					}
				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.orbsOfAtom[soln.atomOfOrb[k]]) {
							val += densityMatrixDeriv.get(l, m) * integralArray[integralcount];
							integralcount++;
						}
					}
				}

				responsematrix.set(j, k, val);
				responsematrix.set(k, j, val);
			}
		}

		SimpleMatrix Roccvirt = soln.CtOcc.mult(responsematrix).mult(soln.CVirt);

		SimpleMatrix Rvec = Roccvirt.elementDivi(soln.Emat);
		Rvec.reshape(NOcc * NVirt, 1);

		return Rvec;
	}

	public static SimpleMatrix computeResponseVectorsThiel(SolutionR soln, SimpleMatrix x) {

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
						sum -= 2 * (soln.Ct.get(i, u) * soln.Ct.get(j + NOcc,
								v) +
								soln.Ct.get(j + NOcc, u) * soln.Ct.get(i, v)) *
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
						element += soln.Ct.get(i, u) * soln.Ct.get(j + NOcc, v) *
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
				p.set(counter, 0, -R.get(counter, 0) + (soln.E.get(j + NOcc) - soln.E.get(i)) * x.get(counter));
				counter++;
			}
		}


		return p;
	}

	private static int numIterable(boolean[] iterable) {
		int count = 0;

		for (boolean value : iterable) {
			if (!value) count++;
		}

		return count;
	}
}
