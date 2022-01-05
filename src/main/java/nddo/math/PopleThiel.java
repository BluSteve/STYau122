package nddo.math;

import nddo.solution.SolutionR;
import nddo.solution.SolutionU;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.SpecializedOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.ejml.dense.row.CommonOps_DDRM.multRows;
import static tools.Utils.mag;

public class PopleThiel { // stop trying to make this faster!!!!!
	public static SimpleMatrix[] toMO(SimpleMatrix CtOcc, SimpleMatrix CVirt, SimpleMatrix[] fockderivstatic) {
		SimpleMatrix[] res = new SimpleMatrix[fockderivstatic.length];

		for (int i = 0; i < res.length; i++) {
			res[i] = CtOcc.mult(fockderivstatic[i]).mult(CVirt);
		}

		return res;
	}

	public static SimpleMatrix[] pople(SolutionR soln, SimpleMatrix[] fockderivstatic) {
		// Pople alg will solve any equation of the form (1-D)x=B
		int NOcc = soln.rm.nOccAlpha;
		int NVirt = soln.rm.nVirtAlpha;
		int length = fockderivstatic.length;
		int nonv = soln.rm.nonvAlpha;

		// array initialization
		SimpleMatrix[] xarray = new SimpleMatrix[length];
		SimpleMatrix[] barray = new SimpleMatrix[length]; // B tilde
		SimpleMatrix[] parray = new SimpleMatrix[length]; // trial R / (ej - ei) as well as D * B tilde
		SimpleMatrix[] Farray = new SimpleMatrix[length]; // F in MO basis divided by (ej - ei)
		SimpleMatrix[] rarray = new SimpleMatrix[length]; // x - F
		double[] oldrMags = new double[length];
		Arrays.fill(oldrMags, 1);

		for (int a = 0; a < length; a++) {
			rarray[a] = new SimpleMatrix(nonv, 1);
		}

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
			SimpleMatrix f = fockderivstatic[a].elementDiv(soln.Emat); // divided by (ej - ei)
			f.reshape(nonv, 1);

			multRows(Darr, f.getDDRM());
			barray[a] = f;
			Farray[a] = f.copy();
			F.setColumn(a, 0, barray[a].getDDRM().data);
		}

		// main loop
		boolean[] iterable = new boolean[length];
		boolean[] looselyIterable = new boolean[length];

		// 0: B, 2: B-P
		List<SimpleMatrix[]> prevs = new ArrayList<>();
		List<Double> dots = new ArrayList<>();

		// responseMatrix remains the same size so pre-initialized
		SimpleMatrix responseMatrix = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		SimpleMatrix alpha = null;
		int prevSize = 0;

		bigLoop:
		while (numIterable(iterable) > 0) {
			// orthogonalize barray
			for (int i = 1; i < barray.length; i++) {
				for (int j = 0; j < i; j++) {
					barray[i].plusi(-barray[i].dot(barray[j]) / barray[j].dot(barray[j]), barray[j]);
				}
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix b = barray[i].copy();
				multRows(Dinvarr, b.getDDRM());
				SimpleMatrix p = computeResponseVectorsPople(soln, b, responseMatrix); // P = D * B
				multRows(Darr, p.getDDRM());
				parray[i] = p;

				SimpleMatrix[] prev = new SimpleMatrix[2];
				prev[0] = barray[i]; // original barray object here
				dots.add(barray[i].dot(barray[i]));

				prev[1] = barray[i].minus(parray[i]);

				prevs.add(prev);
			}

			prevSize = prevs.size();

			for (int i = 0; i < length; i++) {
				SimpleMatrix newb = parray[i].copy();
				// orthogonalize against all previous Bs
				for (int j = 0; j < prevSize; j++) {
					SimpleMatrix[] prev = prevs.get(j);
					double num = prev[0].dot(parray[i]) / dots.get(j);

					newb.plusi(-num, prev[0]);
				}

				barray[i] = newb; // new barray object created
			}

			SimpleMatrix Bt = new SimpleMatrix(prevSize, nonv);
			SimpleMatrix BminusP = new SimpleMatrix(nonv, prevSize);

			for (int i = 0; i < prevSize; i++) {
				Bt.setRow(i, 0, prevs.get(i)[0].getDDRM().data);
				BminusP.setColumn(i, 0, prevs.get(i)[1].getDDRM().data);
			}

			SimpleMatrix rhs = Bt.mult(F);
			SimpleMatrix lhs = Bt.mult(BminusP);
			// alpha dimensions are prevBs x length
			alpha = lhs.solve(rhs);

			// reset r array
			for (int a = 0; a < length; a++) {
				rarray[a].zero();
			}

			for (int i = 0; i < prevSize; i++) {
				for (int j = 0; j < length; j++) {
					rarray[j].plusi(alpha.get(i, j), prevs.get(i)[1]);
				}
			}

			for (int j = 0; j < length; j++) {
				double mag = SpecializedOps_DDRM.diffNormF_fast(rarray[j].getDDRM(), Farray[j].getDDRM());

				if (mag > oldrMags[j] || mag != mag) { // unstable
					if (numIterable(looselyIterable) == 0) {
						soln.rm.getLogger().warn("Slight numerical instability detected; " +
								"returning lower precision values. rarray mag = {}", mag);
						break bigLoop; // i.e. mag is tolerable
					}

					soln.rm.getLogger().warn("Pople algorithm fails; reverting to Thiel algorithm...");
					return thiel(soln, fockderivstatic);
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

	public static SimpleMatrix[] thiel(SolutionR soln, SimpleMatrix[] fockderivstatic) {
		int NOcc = soln.rm.nOccAlpha;
		int NVirt = soln.rm.nVirtAlpha;
		int length = fockderivstatic.length;
		int nonv = soln.rm.nonvAlpha;

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
			SimpleMatrix f = fockderivstatic[a].copy(); // convert AO to MO basis
			f.reshape(nonv, 1);
			multRows(Darr, f.getDDRM());

			xarray[a] = new SimpleMatrix(nonv, 1);
			rarray[a] = f;
			dirs[a] = f.copy();
		}

		SimpleMatrix responseMatrix = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		ArrayList<SimpleMatrix> d = new ArrayList<>();
		ArrayList<SimpleMatrix> p = new ArrayList<>();

		int numIt = 0;
		while (Utils.numNotNull(rarray) > 0) {
			d.clear();
			p.clear();

			for (int i = 0; i < length; i++) {
				if (rarray[i] != null) {
					d.add(dirs[i]);

					SimpleMatrix sm = computeResponseVectorsThiel(soln, dirs[i], responseMatrix);
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
			} catch (SingularMatrixException ignored) {
				alpha = SimpleMatrix.ones(size, length);
			}

			for (int a = 0; a < length; a++) {
				if (rarray[a] != null) {
					for (int i = 0; i < size; i++) {
						double v = alpha.get(i, a);
						xarray[a].plusi(v, d.get(i));
						rarray[a].plusi(-v, p.get(i));
					}

					if (mag(rarray[a]) < 1E-10) {
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
			} catch (SingularMatrixException ignored) {
				beta = SimpleMatrix.ones(size, length);
			}

			for (int a = 0; a < length; a++) {
				if (rarray[a] != null) {
					dirs[a].setTo(rarray[a]);

					for (int i = 0; i < size; i++) {
						dirs[a].plusi(beta.get(i, a), d.get(i));
					}
				}
			}

			if (numIt++ > 10000) {
				IllegalStateException e = new IllegalStateException("Thiel has failed!");
				soln.rm.getLogger().warn(e);
				throw e;
			}
		}

		return xarray;
	}

	public static SimpleMatrix[] thiel(SolutionU soln, SimpleMatrix[] fockderivstaticalpha,
									   SimpleMatrix[] fockderivstaticbeta) {
		int NOccAlpha = soln.rm.nOccAlpha;
		int NOccBeta = soln.rm.nOccBeta;
		int NVirtAlpha = soln.rm.nVirtAlpha;
		int NVirtBeta = soln.rm.nVirtBeta;
		int length = fockderivstaticalpha.length;
		int nonvAlpha = soln.rm.nonvAlpha;
		int nonvBeta = soln.rm.nonvBeta;
		int nonv = nonvAlpha + nonvBeta;

		SimpleMatrix[] xarray = new SimpleMatrix[length];
		SimpleMatrix[] rarray = new SimpleMatrix[length];
		SimpleMatrix[] dirs = new SimpleMatrix[length];

		double[] Darr = new double[nonv];

		int counter = 0;

		for (int i = 0; i < NOccAlpha; i++) {
			for (int j = 0; j < NVirtAlpha; j++) {
				double e = -soln.Ea.get(i) + soln.Ea.get(NOccAlpha + j);

				Darr[counter] = 1 / Math.sqrt(e);

				counter++;
			}
		}

		for (int i = 0; i < NOccBeta; i++) {
			for (int j = 0; j < NVirtBeta; j++) {
				double e = -soln.Eb.get(i) + soln.Eb.get(NOccBeta + j);

				Darr[counter] = 1 / Math.sqrt(e);

				counter++;
			}
		}

		for (int a = 0; a < xarray.length; a++) {
			SimpleMatrix fa = fockderivstaticalpha[a].copy();
			fa.reshape(nonvAlpha, 1);

			SimpleMatrix fb = fockderivstaticbeta[a].copy();
			fb.reshape(nonvBeta, 1);

			SimpleMatrix f = fa.concatRows(fb);
			multRows(Darr, f.getDDRM());

			xarray[a] = new SimpleMatrix(nonv, 1);
			rarray[a] = f;
			dirs[a] = f.copy();
		}

		SimpleMatrix xmata = new SimpleMatrix(soln.rm.nonvAlpha, 1);
		SimpleMatrix xmatb = new SimpleMatrix(soln.rm.nonvBeta, 1);
		SimpleMatrix Jderiv = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		SimpleMatrix Kaderiv = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		SimpleMatrix Kbderiv = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		ArrayList<SimpleMatrix> d = new ArrayList<>();
		ArrayList<SimpleMatrix> p = new ArrayList<>();

		int numIt = 0;
		while (Utils.numNotNull(rarray) > 0) {
			d.clear();
			p.clear();

			for (int i = 0; i < length; i++) {
				if (rarray[i] != null) {
					d.add(dirs[i]);

					SimpleMatrix sm =
							computeResponseVectorsThiel(soln, dirs[i], xmata, xmatb, Jderiv, Kaderiv, Kbderiv);
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
			} catch (SingularMatrixException ignored) {
				alpha = SimpleMatrix.ones(size, length);
			}

			for (int a = 0; a < length; a++) {
				if (rarray[a] != null) {
					for (int i = 0; i < size; i++) {
						double v = alpha.get(i, a);
						xarray[a].plusi(v, d.get(i));
						rarray[a].plusi(-v, p.get(i));
					}

					if (mag(rarray[a]) < 1E-10) {
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
			} catch (SingularMatrixException ignored) {
				beta = SimpleMatrix.ones(size, length);
			}

			for (int a = 0; a < length; a++) {
				if (rarray[a] != null) {
					dirs[a].setTo(rarray[a]);

					for (int i = 0; i < size; i++) {
						dirs[a].plusi(beta.get(i, a), d.get(i));
					}
				}
			}

			if (numIt++ > 10000) {
				IllegalStateException e = new IllegalStateException("Thiel has failed!");
				soln.rm.getLogger().warn(e);
				throw e;
			}
		}

		return xarray;
	}

	public static SimpleMatrix densityDeriv(SolutionR soln, SimpleMatrix x) {
		x.reshape(soln.rm.nOccAlpha, soln.rm.nVirtAlpha);

		SimpleMatrix mult = soln.COcc.mult(x).mult(soln.CtVirt);

		x.reshape(soln.rm.nonvAlpha, 1);

		return Utils.plusTrans(mult).scalei(-2);
	}

	public static SimpleMatrix[] densityDeriv(SolutionU soln, SimpleMatrix x) {
		SimpleMatrix xmata = x.extractMatrix(0, soln.rm.nonvAlpha, 0, 1);
		xmata.reshape(soln.rm.nOccAlpha, soln.rm.nVirtAlpha);

		SimpleMatrix mult = soln.CaOcc.mult(xmata).mult(soln.CtaVirt);
		SimpleMatrix densityDerivAlpha = Utils.plusTrans(mult).negativei();

		SimpleMatrix xmatb = x.extractMatrix(soln.rm.nonvAlpha, x.numRows(), 0, 1);
		xmatb.reshape(soln.rm.nOccBeta, soln.rm.nVirtBeta);

		mult = soln.CbOcc.mult(xmatb).mult(soln.CtbVirt);
		SimpleMatrix densityDerivBeta = Utils.plusTrans(mult).negativei();

		return new SimpleMatrix[]{densityDerivAlpha, densityDerivBeta};
	}

	public static SimpleMatrix responseMatrix(SolutionR soln, SimpleMatrix densityDeriv) {
		SimpleMatrix responsematrix = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		responseMatrix(soln, densityDeriv, responsematrix);

		return responsematrix;
	}

	public static SimpleMatrix[] responseMatrix(SolutionU soln, SimpleMatrix[] densityDeriv) {
		SimpleMatrix Jderiv = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		SimpleMatrix Kaderiv = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);
		SimpleMatrix Kbderiv = new SimpleMatrix(soln.nOrbitals, soln.nOrbitals);

		return responseMatrix(soln, densityDeriv[0], densityDeriv[1], Jderiv, Kaderiv, Kbderiv);
	}

	private static SimpleMatrix computeResponseVectorsPople(SolutionR soln, SimpleMatrix x,
															SimpleMatrix responseMatrix) {
		// x is B tilde, i.e. a guess
		// p = d * b. this evaluates d * b without finding out what d is, as the latter is slow.

		SimpleMatrix densityDeriv = densityDeriv(soln, x);

		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double val = 0;
				if (j == k) {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						val += densityDeriv.get(l, l) * integralArray[integralcount];
						integralcount++;
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += densityDeriv.get(l, m) * integralArray[integralcount];
								integralcount++;
							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += densityDeriv.get(j, k) * integralArray[integralcount];
					integralcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += densityDeriv.get(l, m) * integralArray[integralcount];
								integralcount++;
							}
						}
					}
				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.orbsOfAtom[soln.atomOfOrb[k]]) {
							val += densityDeriv.get(l, m) * integralArray[integralcount];
							integralcount++;
						}
					}
				}

				responseMatrix.set(j, k, val);
				responseMatrix.set(k, j, val);
			}
		}

		SimpleMatrix Roccvirt = soln.CtOcc.mult(responseMatrix).mult(soln.CVirt);

		SimpleMatrix Rvec = Roccvirt.elementDivi(soln.Emat);
		Rvec.reshape(soln.rm.nOccAlpha * soln.rm.nVirtAlpha, 1);

		return Rvec;
	}

	private static SimpleMatrix computeResponseVectorsThiel(SolutionR soln, SimpleMatrix x,
															SimpleMatrix responseMatrix) {
		x.reshape(soln.rm.nOccAlpha, soln.rm.nVirtAlpha);

		SimpleMatrix mult = soln.COcc.mult(x).mult(soln.CtVirt);
		SimpleMatrix densityDeriv = Utils.plusTrans(mult).scalei(-2);

		responseMatrix(soln, densityDeriv, responseMatrix);

		SimpleMatrix p = soln.Emat.elementMult(x).minusi(soln.CtOcc.mult(responseMatrix).mult(soln.CVirt));
		p.reshape(soln.rm.nonvAlpha, 1);

		x.reshape(soln.rm.nonvAlpha, 1);

		return p;
	}

	private static SimpleMatrix computeResponseVectorsThiel(SolutionU soln, SimpleMatrix x, SimpleMatrix xmata,
															SimpleMatrix xmatb, SimpleMatrix Jderiv,
															SimpleMatrix Kaderiv, SimpleMatrix Kbderiv) {
		CommonOps_DDRM.extract(x.getDDRM(), 0, soln.rm.nonvAlpha, 0, 1, xmata.getDDRM());
		xmata.reshape(soln.rm.nOccAlpha, soln.rm.nVirtAlpha);

		SimpleMatrix mult = soln.CaOcc.mult(xmata).mult(soln.CtaVirt);
		SimpleMatrix densityDerivAlpha = Utils.plusTrans(mult).negativei();

		CommonOps_DDRM.extract(x.getDDRM(), soln.rm.nonvAlpha, x.numRows(), 0, 1, xmatb.getDDRM());
		xmatb.reshape(soln.rm.nOccBeta, soln.rm.nVirtBeta);

		mult = soln.CbOcc.mult(xmatb).mult(soln.CtbVirt);
		SimpleMatrix densityDerivBeta = Utils.plusTrans(mult).negativei();

		SimpleMatrix[] sms = responseMatrix(soln, densityDerivAlpha, densityDerivBeta, Jderiv, Kaderiv, Kbderiv);

		SimpleMatrix pa = soln.Eamat.elementMult(xmata).minusi(soln.CtaOcc.mult(sms[0]).mult(soln.CaVirt));
		pa.reshape(soln.rm.nonvAlpha, 1);
		SimpleMatrix pb = soln.Ebmat.elementMult(xmatb).minusi(soln.CtbOcc.mult(sms[1]).mult(soln.CbVirt));
		pb.reshape(soln.rm.nonvBeta, 1);

		xmata.reshape(soln.rm.nonvAlpha, 1);
		xmatb.reshape(soln.rm.nonvBeta, 1);

		return pa.concatRows(pb);
	}

	private static SimpleMatrix responseMatrix(SolutionR soln, SimpleMatrix densityDeriv, SimpleMatrix target) {
		double[] integralArray = soln.integralArray;

		int integralcount = 0;

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double val = 0;
				if (j == k) {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						val += densityDeriv.get(l, l) * integralArray[integralcount];
						integralcount++;
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += densityDeriv.get(l, m) * integralArray[integralcount];
								integralcount++;
							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += densityDeriv.get(j, k) * integralArray[integralcount];
					integralcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += densityDeriv.get(l, m) * integralArray[integralcount];
								integralcount++;
							}
						}
					}
				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.orbsOfAtom[soln.atomOfOrb[k]]) {
							val += densityDeriv.get(l, m) * integralArray[integralcount];
							integralcount++;
						}
					}
				}

				target.set(j, k, val);
				target.set(k, j, val);
			}
		}

		return target;
	}

	private static SimpleMatrix[] responseMatrix(SolutionU soln, SimpleMatrix densityDerivAlpha,
												 SimpleMatrix densityDerivBeta, SimpleMatrix Jderiv,
												 SimpleMatrix Kaderiv, SimpleMatrix Kbderiv) {
		int Jcount = 0;
		int Kcount = 0;

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double val = 0;
				if (j == k) {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						val += (densityDerivAlpha.get(l, l) + densityDerivBeta.get(l, l)) *
								soln.integralArrayCoulomb[Jcount];
						Jcount++;
					}

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += (densityDerivAlpha.get(l, m) + densityDerivBeta.get(l, m)) *
										soln.integralArrayCoulomb[Jcount];
								Jcount++;
							}
						}
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					val += (densityDerivAlpha.get(j, k) + densityDerivBeta.get(j, k)) *
							soln.integralArrayCoulomb[Jcount];
					Jcount++;

					for (int l : soln.missingOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.missingOfAtom[soln.atomOfOrb[j]]) {
							if (soln.atomOfOrb[l] == soln.atomOfOrb[m]) {
								val += (densityDerivAlpha.get(l, m) + densityDerivBeta.get(l, m)) *
										soln.integralArrayCoulomb[Jcount];
								Jcount++;
							}
						}
					}
				}

				Jderiv.set(j, k, val);
				Jderiv.set(k, j, val);
			}
		}

		for (int j = 0; j < soln.nOrbitals; j++) {
			for (int k = j; k < soln.nOrbitals; k++) {
				double vala = 0;
				double valb = 0;
				if (j == k) {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						vala += densityDerivAlpha.get(l, l) * soln.integralArrayExchange[Kcount];
						valb += densityDerivBeta.get(l, l) * soln.integralArrayExchange[Kcount];
						Kcount++;
					}
				}
				else if (soln.atomOfOrb[j] == soln.atomOfOrb[k]) {
					vala += densityDerivAlpha.get(j, k) * soln.integralArrayExchange[Kcount];
					valb += densityDerivBeta.get(j, k) * soln.integralArrayExchange[Kcount];
					Kcount++;

				}
				else {
					for (int l : soln.orbsOfAtom[soln.atomOfOrb[j]]) {
						for (int m : soln.orbsOfAtom[soln.atomOfOrb[k]]) {
							vala += densityDerivAlpha.get(l, m) * soln.integralArrayExchange[Kcount];
							valb += densityDerivBeta.get(l, m) * soln.integralArrayExchange[Kcount];
							Kcount++;
						}
					}
				}

				Kaderiv.set(j, k, vala);
				Kaderiv.set(k, j, vala);

				Kbderiv.set(j, k, valb);
				Kbderiv.set(k, j, valb);
			}
		}

		return new SimpleMatrix[]{Kaderiv.plusi(Jderiv), Kbderiv.plusi(Jderiv)};
	}

	private static int numIterable(boolean[] iterable) {
		int count = 0;

		for (boolean value : iterable) {
			if (!value) count++;
		}

		return count;
	}
}
