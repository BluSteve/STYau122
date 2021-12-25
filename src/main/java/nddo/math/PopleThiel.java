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

			for (int i = 0; i < length; i++) {
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
		double[] oldrMags = new double[length];
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
		int prevSize = 0;
		SimpleMatrix rhs = null;
		SimpleMatrix lhs = null;
		SimpleMatrix Bt2 = new SimpleMatrix(length, nonv); // last 15
		SimpleMatrix P2 = new SimpleMatrix(nonv, length);
		int n = 1;

		bigLoop:
		while (numIterable(iterable) > 0) {
			// orthogonalize barray
			for (int i = 1; i < barray.length; i++) {
				for (int j = 0; j < i; j++) {
					barray[i].plusi(-barray[i].dot(barray[j]) / barray[j].dot(barray[j]), barray[j]);
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

			prevSize = prevs.size();

			for (int i = 0; i < length; i++) {
				SimpleMatrix newb = parray[i].copy();

				// orthogonalize against all previous Bs
				for (int j = 0; j < prevSize; j++) {
					SimpleMatrix[] prev = prevs.get(j);
					SimpleMatrix transpose = prev[1];
					double num = transpose.mult(parray[i]).get(0) / dots.get(j);

					newb.plusi(num, prev[2]);
				}

				barray[i] = newb; // new barray object created
			}

			// convert prevBs and prevPs into matrix form, transposed
			int prevL = (n - 1) * length;

			// everything but last 15
			SimpleMatrix Bt1 = new SimpleMatrix(prevL, nonv);
			SimpleMatrix P = new SimpleMatrix(nonv, prevs.size());

			for (int i = 0; i < prevs.size(); i++) {
				if (i >= prevL) {
					Bt2.setRow(i - prevL, 0, prevs.get(i)[0].getDDRM().data);
					P2.setColumn(i - prevL, 0, prevs.get(i)[4].getDDRM().data);
				}
				else {
					Bt1.setRow(i, 0, prevs.get(i)[0].getDDRM().data);
				}
				P.setColumn(i, 0, prevs.get(i)[4].getDDRM().data);
			}

			if (rhs == null) rhs = Bt2.mult(F);
			else rhs = rhs.combine(prevL, 0, Bt2.mult(F));

			if (lhs == null) lhs = Bt2.mult(P);
			else {
				SimpleMatrix topright = Bt1.mult(P2);
				SimpleMatrix bottom = Bt2.mult(P);

				int nl = n * length;
				SimpleMatrix newlhs = new SimpleMatrix(nl, nl);
				newlhs.insertIntoThis(0, 0, lhs);
				newlhs.insertIntoThis(0, prevL, topright);
				newlhs.insertIntoThis(prevL, 0, bottom);
				lhs = newlhs;
			}

			alpha = lhs.solve(rhs);

			// reset r array
			for (int a = 0; a < length; a++) {
				rarray[a] = new SimpleMatrix(nonv, 1);
			}

			for (int i = 0; i < prevSize; i++) {
				for (int j = 0; j < length; j++) {
					rarray[j].plusi(alpha.get(i, j), prevs.get(i)[4]);
				}
			}

			for (int j = 0; j < length; j++) {
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
			n++;
		}

		for (int i = 0; i < length; i++) {
			xarray[i] = new SimpleMatrix(nonv, 1);
		}

		for (int i = 0; i < prevSize; i++) {
			for (int j = 0; j < length; j++) {
				xarray[j].plusi(alpha.get(i, j), prevs.get(i)[0]);
			}
		}

		for (int j = 0; j < length; j++) {
			multRows(Dinvarr, xarray[j].getDDRM());
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

		SimpleMatrix mult = soln.COcc.mult(xmat).mult(soln.CtVirt);
		SimpleMatrix densityMatrixDeriv = mult.plusi(mult.transpose()).scalei(-2);

		SimpleMatrix responsematrix = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
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

		for (int a = 0; a < length; a++) {
			SimpleMatrix f = soln.CtOcc.mult(fockderivstatic[a]).mult(soln.CVirt); // convert AO to MO basis
			f.reshape(nonv, 1);
			multRows(Darr, f.getDDRM());

			xarray[a] = new SimpleMatrix(nonv, 1);
			rarray[a] = f;
			dirs[a] = f;
		}

		ArrayList<SimpleMatrix> d = new ArrayList<>();
		ArrayList<SimpleMatrix> p = new ArrayList<>();

		while (Utils.numNotNull(rarray) > 0) {
			d.clear();
			p.clear();

			for (int i = 0; i < length; i++) {
				if (rarray[i] != null) {
					d.add(dirs[i]);

					SimpleMatrix sm = computeResponseVectorsThiel(soln, dirs[i]);
					multRows(Darr, sm.getDDRM());
					p.add(sm);
				}
			}

			int size = p.size(); // actually iterable size

			SimpleMatrix solver = new SimpleMatrix(size, size);
			SimpleMatrix rhsvec = new SimpleMatrix(size, length);
			double[] arrrhs = new double[size];

			for (int a = 0; a < length; a++) {
				if (rarray[a] != null) {
					for (int i = 0; i < size; i++) {
						arrrhs[i] = 2 * rarray[a].dot(d.get(i));
					}

					rhsvec.setColumn(a, 0, arrrhs);
				}
			}

			for (int i = 0; i < size; i++) {
				for (int j = i; j < size; j++) {
					double val2 = p.get(j).dot(d.get(i)) + p.get(i).dot(d.get(j));

					solver.set(i, j, val2);
					solver.set(j, i, val2);
				}
			}

			SimpleMatrix alpha;
			try {
				alpha = solver.solve(rhsvec);
			} catch (SingularMatrixException e) {
				alpha = SimpleMatrix.ones(size, length);
			}

			for (int a = 0; a < length; a++) {
				if (rarray[a] != null) {
					for (int i = 0; i < size; i++) {
						double v = alpha.get(i, a);
						xarray[a].plusi(v, d.get(i));
						rarray[a].plusi(-v, p.get(i));
					}

					if (mag(rarray[a]) < 1E-6) { // todo change this if you want
						rarray[a] = null;
					}
				}
			}

			solver = new SimpleMatrix(size, size);

			for (int a = 0; a < length; a++) {
				if (rarray[a] != null) {
					for (int i = 0; i < size; i++) {
						arrrhs[i] = -rarray[a].dot(p.get(i));
					}

					rhsvec.setColumn(a, 0, arrrhs);
				}
			}

			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					solver.set(i, j, d.get(j).dot(p.get(i)));
				}
			}

			SimpleMatrix beta;
			try {
				beta = solver.solve(rhsvec);
			} catch (SingularMatrixException e) {
				beta = SimpleMatrix.ones(size, length);
			}

			for (int a = 0; a < length; a++) {
				if (rarray[a] != null) {
					dirs[a] = rarray[a].copy();

					for (int i = 0; i < size; i++) {
						dirs[a].plusi(beta.get(i, a), d.get(i));
					}
				}
			}
		}

		return xarray;
	}

	public static SimpleMatrix computeResponseVectorsThiel(SolutionR soln, SimpleMatrix x) {
		int NOcc = soln.rm.nOccAlpha;
		int NVirt = soln.rm.nVirtAlpha;

		SimpleMatrix xmat = x.copy();
		xmat.reshape(NOcc, NVirt);

		SimpleMatrix mult = soln.COcc.mult(xmat).mult(soln.CtVirt);
		SimpleMatrix densityMatrixDeriv = mult.plusi(mult.transpose()).scalei(-2);

		SimpleMatrix responsematrix = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
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

		SimpleMatrix p = soln.Emat.elementMult(xmat).minusi(soln.CtOcc.mult(responsematrix).mult(soln.CVirt));
		p.reshape(NOcc * NVirt, 1);

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