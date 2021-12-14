package testing;

import nddo.solution.SolutionR;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import frontend.json.InputHandler;
import tools.Utils;

import java.util.ArrayList;
import java.util.List;

import static nddo.geometry.GeometrySecondDerivative.computeResponseVectorsPople;

public class Testing {
	public static void main(String[] args) throws Exception {
//		RawInput ri = InputHandler.processInput("subset");
//		RunnableMolecule rm = ri.molecules[0];
//
//		NDDOParams[] npMap = Utils.getNpMap(ri);
//		Solution s = Solution.of(rm, State.getConverter().convert(rm.atoms));
//
//		System.out.println("s = " + s);
		InputHandler.txtToText();
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
}
