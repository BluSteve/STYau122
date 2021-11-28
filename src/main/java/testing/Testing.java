package testing;

import nddo.geometry.GeometryDerivative;
import nddo.mndo.MNDOAtom;
import nddo.mndo.MNDOParams;
import nddo.solution.SolutionR;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import runcycle.input.RawMolecule;
import scf.AtomHandler;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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

		int initc = 2 * (int) Math.ceil(length * Math.log(length));
		// 0: B, 1: Bt, 2: Bn, 3: P, 4: BmP
		List<SimpleMatrix[]> prevs = new ArrayList<>(initc);

		SimpleMatrix rhs = null;
		SimpleMatrix lhs = null;
		SimpleMatrix Bt2 = new SimpleMatrix(length, nonv); // last 15
		SimpleMatrix P2 = new SimpleMatrix(nonv, length);
		int n = 1;

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
				for (SimpleMatrix[] prev : prevs) {
					SimpleMatrix prevB = prev[0];
					SimpleMatrix transpose = prev[1];
					double num = transpose.mult(parray[i]).get(0) /
							prevB.dot(prevB);

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

			SimpleMatrix topright = Bt1.mult(P2);
			SimpleMatrix bottom = Bt2.mult(P);

			if (rhs == null) rhs = Bt2.mult(F);
			else rhs = rhs.combine(prevL, 0, Bt2.mult(F));

			if (lhs == null) lhs = Bt2.mult(P);
			else {
				int nl = n * length;
				SimpleMatrix newlhs = new SimpleMatrix(nl, nl);
				newlhs.insertIntoThis(0, 0, lhs);
				newlhs.insertIntoThis(0, prevL, topright);
				newlhs.insertIntoThis(prevL, 0, bottom);
				lhs = newlhs;
			}

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
		System.out.println(Arrays.toString(getxarrayPople(s, matrices[1])));

		System.exit(0);
	}
}
