package testing;

import nddo.geometry.GeometryDerivative;
import nddo.mndo.MNDOAtom;
import nddo.mndo.MNDOParams;
import nddo.solution.SolutionR;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import runcycle.input.RawMolecule;
import scf.AtomHandler;
import tools.Utils;

import java.util.ArrayList;

import static nddo.geometry.GeometrySecondDerivative.computeResponseVectorsPople;

public class Testing {
	public static void main(String[] args) throws Exception {
		testMain();
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
		for (int i = 0; i < length; i++) {
			parray[i] = new SimpleMatrix(nonv, 1);
		}
		SimpleMatrix[] Farray = new SimpleMatrix[length];
		SimpleMatrix[] rarray = new SimpleMatrix[length];

		// configure preconditioners
		SimpleMatrix D = new SimpleMatrix(nonv, nonv, DMatrixSparseCSC.class);
		SimpleMatrix Dinv =
				new SimpleMatrix(nonv, nonv, DMatrixSparseCSC.class);

		int counter = 0;
		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = -soln.E.get(i) + soln.E.get(NOcc + j);

				D.set(counter, counter, Math.pow(e, -0.5));
				Dinv.set(counter, counter, Math.pow(e, 0.5));
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

			barray[a] = D.mult(f);
			Farray[a] = barray[a].copy();
			F.setColumn(a, 0, barray[a].getDDRM().data);
		}

		// main loop
		int[] iterable = new int[length];
		ArrayList<SimpleMatrix> prevBs = new ArrayList<>();
		ArrayList<SimpleMatrix> prevBTs = new ArrayList<>();
		ArrayList<SimpleMatrix> prevPs = new ArrayList<>();
		ArrayList<SimpleMatrix> prevBmP = new ArrayList<>();
		ArrayList<SimpleMatrix> prevBn = new ArrayList<>();
		SimpleMatrix rhs2 = null;
		SimpleMatrix lhs2 = null;
		int n = 1;

		while (Utils.numIterable(iterable) > 0) {
			Utils.orthogonalise(barray);

			for (int i = 0; i < length; i++) {
				prevBs.add(barray[i]); // original barray object here
				prevBTs.add(barray[i].transpose());
				prevBn.add(barray[i].negative());

				// parray[i] stays the same object throughout
				CommonOps_DDRM.mult(D.getDDRM(),
						computeResponseVectorsPople
								(Dinv.mult(barray[i]), soln).getDDRM(),
						parray[i].getDDRM());

				prevPs.add(parray[i].copy());
				prevBmP.add(barray[i].minus(parray[i]));
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix newb = parray[i].copy();

				// orthogonalize against all previous Bs
				for (int j = 0; j < prevBs.size(); j++) {
					SimpleMatrix prevB = prevBs.get(j);
					SimpleMatrix transpose = prevBTs.get(j);
					double num = transpose.mult(parray[i]).get(0) /
							prevB.dot(prevB);

					newb.plusi(num, prevBn.get(j));
				}

				barray[i] = newb; // new barray object created
			}

			// convert prevBs and prevPs into matrix form, transposed
			int prevL = (n - 1) * length;

			// everything but last 15
			SimpleMatrix Bt1 = new SimpleMatrix(prevL, nonv);
			SimpleMatrix Bt2 = new SimpleMatrix(length, nonv); // last 15

			SimpleMatrix P = new SimpleMatrix(nonv, prevPs.size());
			SimpleMatrix P2 = new SimpleMatrix(nonv, length);

			for (int i = 0; i < prevBs.size(); i++) {
				if (i >= prevL) {
					Bt2.setRow(i - prevL, 0, prevBs.get(i).getDDRM().data);
					P2.setColumn(i - prevL, 0, prevBmP.get(i).getDDRM().data);
				}
				else {
					Bt1.setRow(i, 0, prevBs.get(i).getDDRM().data);
				}
				P.setColumn(i, 0, prevBmP.get(i).getDDRM().data);
			}

			SimpleMatrix topright = Bt1.mult(P2);
			SimpleMatrix bottom = Bt2.mult(P);

			if (rhs2 == null) rhs2 = Bt2.mult(F);
			else rhs2 = rhs2.combine(prevL, 0, Bt2.mult(F));

			if (lhs2 == null) lhs2 = Bt2.mult(P);
			else {
				SimpleMatrix newlhs = new SimpleMatrix(n * length, n*length);
				newlhs.insertIntoThis(0, 0, lhs2);
				newlhs.insertIntoThis(0, prevL, topright);
				newlhs.insertIntoThis(prevL, 0, bottom);
				lhs2 = newlhs;
			}

			// alpha dimensions are prevBs x length
			SimpleMatrix alpha = lhs2.solve(rhs2);

			// reset r and x array
			for (int a = 0; a < length; a++) {
				rarray[a] = new SimpleMatrix(nonv, 1);
				xarray[a] = new SimpleMatrix(nonv, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {
					// B with tilde
					rarray[j].plusi(alpha.get(i, j), prevBmP.get(i));
					xarray[j].plusi(alpha.get(i, j), prevBs.get(i));
				}
			}

			for (int j = 0; j < alpha.numCols(); j++) {
				// B0 is Farray, no tilde
				rarray[j].minusi(Farray[j]);
				xarray[j] = Dinv.mult(xarray[j]);

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

			n++;
		}

		return xarray;
	}

	public static void testMain() throws InterruptedException {
		AtomHandler.populateAtoms();
		MNDOParams h = new MNDOParams(2.92397599125172,
				-6.222578482830868, 0.0, -12.200235077462583, 0.0,
				1.0693232546199132, 0.0, -13.00142320543855, 12.848,
				0.0, 0.0, 0.0, 0.0);
		MNDOParams c = new MNDOParams(2.5572499654157435, -18.854021376560777,
				-8.377666892780198, -52.57072065877964, -39.05266019981942,
				1.838438013363027, 1.805140784089995, -120.60738371097112,
				12.23, 11.47,
				2.43, 11.08, 9.84);
		MNDOAtom atom1 = new MNDOAtom(AtomHandler.atomsMap.get("H"),
				new double[]{0.635 * Utils.bohr, 0.639 * Utils.bohr,
						0.635 * Utils.bohr},
				h);
		MNDOAtom atom2 = new MNDOAtom(AtomHandler.atomsMap.get("H"),
				new double[]{0.635 * Utils.bohr, -0.635 * Utils.bohr,
						-0.639 * Utils.bohr}, h);
		MNDOAtom atom3 = new MNDOAtom(AtomHandler.atomsMap.get("H"),
				new double[]{-0.639 * Utils.bohr, -0.635 * Utils.bohr,
						0.635 * Utils.bohr}, h);
		MNDOAtom atom4 = new MNDOAtom(AtomHandler.atomsMap.get("H"),
				new double[]{-0.639 * Utils.bohr, 0.639 * Utils.bohr,
						-0.639 * Utils.bohr}, h);
		MNDOAtom carbon = new MNDOAtom(AtomHandler.atomsMap.get("C"),
				new double[]{-0.0021 * Utils.bohr, 0.0021 * Utils.bohr,
						-0.0021 * Utils.bohr}, c);

		MNDOAtom[] atoms = new MNDOAtom[]{atom1, atom2, atom3, atom4, carbon};

		MNDOAtom[] exp = new MNDOAtom[]{
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{0.6304 * Utils.bohr, 0.6304 * Utils.bohr,
								0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{0.6304 * Utils.bohr, -0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{-0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr,
								0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{-0.6304 * Utils.bohr, 0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("C"),
						new double[]{0, 0, 0}, c)};
		MNDOAtom[] exp1 = new MNDOAtom[]{
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{0.6304 * Utils.bohr, 0.6304 * Utils.bohr,
								0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{0.6304 * Utils.bohr, -0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{-0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr,
								0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("H"),
						new double[]{-0.6304 * Utils.bohr, 0.6304 * Utils.bohr,
								-0.6304 * Utils.bohr}, h),
				new MNDOAtom(AtomHandler.atomsMap.get("C"),
						new double[]{0, 0, 0}, c)};
		double[] datum = new double[]{-17.9, 0, 13.6};

		RawMolecule rm = new RawMolecule();
		rm.nIntegrals = 212;
		rm.nElectrons = 8;
		rm.nOrbitals = 8;
		rm.atomicNumbers = new int[]{1, 1, 1, 1, 6};
		rm.charge = 0;
		rm.mult = 0;
		rm.name = "C1H4";

		SolutionR s = new SolutionR(atoms, rm).compute();
		SimpleMatrix[][] matrices = GeometryDerivative.gradientRoutine(s);
		getxarrayPople(s, matrices[1]);

		System.exit(0);
	}
}
