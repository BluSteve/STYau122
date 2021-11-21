import nddoparam.GeometryDerivative;
import nddoparam.SolutionR;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.simple.SimpleMatrix;
import runcycle.input.RawMolecule;
import scf.AtomHandler;
import tools.NanoStopWatch;
import tools.Utils;

import java.util.ArrayList;
import java.util.concurrent.TimeUnit;

import static nddoparam.GeometrySecondDerivative.*;

public class Testing {
	public static void main(String[] args) {
		try {
			testMain();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static SimpleMatrix[] getxarray(SolutionR soln,
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
		SimpleMatrix D = new SimpleMatrix(nonv, nonv, DMatrixSparseCSC.class);
		SimpleMatrix Dinv = new SimpleMatrix(nonv, nonv, DMatrixSparseCSC.class);

		int counter = 0;
		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = (-soln.E.get(i) + soln.E.get(NOcc + j));

				D.set(counter, counter, Math.pow(e, -0.5));
				Dinv.set(counter, counter, Math.pow(e, 0.5));
				counter++;
			}
		}

		// convert AO to MO basis
		for (int a = 0; a < length; a++) {
			SimpleMatrix F = new SimpleMatrix(nonv, 1);

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

					element = element / (soln.E.get(j + NOcc) - soln.E.get(i));

					F.set(count, 0, element);

					count++;
				}
			}

			F = D.mult(F);

			barray[a] = F;
			Farray[a] = F;
		}

		// main loop
		int[] iterable = new int[length];
		ArrayList<SimpleMatrix> prevBs = new ArrayList<>();
		ArrayList<SimpleMatrix> prevPs = new ArrayList<>();

		// convert Farray into matrix form
		SimpleMatrix F = new SimpleMatrix(nonv, length);
		for (int i = 0; i < Farray.length; i++) {
			F.setColumn(i, 0, Farray[i].getDDRM().data);
		}

		while (Utils.numIterable(iterable) > 0) {
			orthogonalise(barray);

			for (int i = 0; i < length; i++) {
				prevBs.add(barray[i]);

				parray[i] = D.mult(computeResponseVectorsPople(
						Dinv.mult(barray[i]), soln));

				prevPs.add(parray[i]);
			}

			for (int i = 0; i < length; i++) {
				SimpleMatrix newb = parray[i];

				// orthogonalize against all previous Bs
				for (SimpleMatrix prevB : prevBs) {
					double num = prevB.transpose().mult(parray[i]).get(0) /
							prevB.transpose().mult(prevB).get(0);

					newb = newb.minus(prevB.scale(num));
				}

				barray[i] = newb;
			}

			// convert prevBs and prevPs into matrix form, transposed
			SimpleMatrix Bt = new SimpleMatrix(prevBs.size(), nonv);
			SimpleMatrix P = new SimpleMatrix(nonv, prevPs.size());
			for (int i = 0; i < prevBs.size(); i++) {
				Bt.setRow(i, 0, prevBs.get(i).getDDRM().data);
				P.setColumn(i, 0,
						prevBs.get(i).minus(prevPs.get(i)).getDDRM().data);
			}

			SimpleMatrix lhs = Bt.mult(P);
			SimpleMatrix rhs = Bt.mult(F);
			SimpleMatrix alpha = lhs.solve(rhs);

			// reset r and x array
			for (int a = 0; a < length; a++) {
				rarray[a] = new SimpleMatrix(nonv, 1);
				xarray[a] = new SimpleMatrix(nonv, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {
					// B with tilde
					rarray[j] = rarray[j].plus(
							prevBs.get(i)
									.minus(prevPs.get(i))
									.scale(alpha.get(i, j))
					);

					xarray[j] = xarray[j].plus(
							prevBs.get(i)
									.scale(alpha.get(i, j))
					);
				}
			}

			for (int j = 0; j < alpha.numCols(); j++) {
				// B0 is Farray, no tilde
				rarray[j] = rarray[j].minus(Farray[j]);
				xarray[j] = Dinv.mult(xarray[j]);

				double rMag = mag(rarray[j]);
				if (rMag < 1E-7) {
					iterable[j] = 1;
				}
				else if (Double.isNaN(rMag)) {
					System.err.println(
							"Pople algorithm fails; reverting to Thiel " +
									"algorithm (don't panic)...");
					return null;
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

		SimpleMatrix[] xarray = new SimpleMatrix[fockderivstatic.length];
		SimpleMatrix[] barray = new SimpleMatrix[fockderivstatic.length];
		SimpleMatrix[] parray = new SimpleMatrix[fockderivstatic.length];
		SimpleMatrix[] Farray = new SimpleMatrix[fockderivstatic.length];
		SimpleMatrix[] rarray = new SimpleMatrix[fockderivstatic.length];

		double[] arrpreconditioner = new double[NOcc * NVirt];
		double[] arrpreconditionerinv = new double[NOcc * NVirt];

		int counter = 0;

		for (int i = 0; i < NOcc; i++) {
			for (int j = 0; j < NVirt; j++) {
				double e = (-soln.E.get(i) + soln.E.get(NOcc + j));

				arrpreconditioner[counter] = Math.pow(e, -0.5);
				arrpreconditionerinv[counter] = Math.pow(e, 0.5);
				counter++;
			}
		}

		SimpleMatrix D = SimpleMatrix.diag(arrpreconditioner);

		SimpleMatrix Dinv = SimpleMatrix.diag(arrpreconditionerinv);

		for (int a = 0; a < xarray.length; a++) {
			SimpleMatrix F = new SimpleMatrix(NOcc * NVirt, 1);

			int count1 = 0;

			for (int i = 0; i < NOcc; i++) { // kappa
				for (int j = 0; j < NVirt; j++) { // i

					double element = 0;

					for (int u = 0; u < soln.orbitals.length; u++) {
						for (int v = 0; v < soln.orbitals.length; v++) {
							element +=
									soln.C.get(i, u) * soln.C.get(j + NOcc,
											v) *
											fockderivstatic[a].get(u, v);
						}
					}

					element = element / (soln.E.get(j + NOcc) - soln.E.get(i));

					F.set(count1, 0, element);

					count1++;
				}
			}

			F = D.mult(F);

			xarray[a] = new SimpleMatrix(NOcc * NVirt, 1);
			rarray[a] = new SimpleMatrix(NOcc * NVirt, 1);
			barray[a] = F;
			Farray[a] = F;
		}


		if (barray[0].numRows() == 0) {
			SimpleMatrix[] densityderivs =
					new SimpleMatrix[fockderivstatic.length];

			for (int i = 0; i < densityderivs.length; i++) {
				densityderivs[i] = new SimpleMatrix(0, 0);
			}

			return densityderivs;
		}


		ArrayList<SimpleMatrix> prevBs = new ArrayList<>();

		ArrayList<SimpleMatrix> prevPs = new ArrayList<>();

		int[] iterable = new int[barray.length];

		SimpleMatrix F =
				new SimpleMatrix(NOcc * NVirt, Farray.length);

		for (int i = 0; i < Farray.length; i++) {
			F.setColumn(i, 0, Farray[i].getDDRM().data);
		}

		while (Utils.numIterable(iterable) > 0) {
			for (int number = 0; number < 1; number++) {
				orthogonalise(barray);


				for (int i = 0; i < barray.length; i++) {
					prevBs.add(barray[i]);

					parray[i] = D.mult(computeResponseVectorsPople(
							Dinv.mult(barray[i]),
							soln));

					prevPs.add(parray[i]);
				}

				for (int i = 0; i < barray.length; i++) {
					SimpleMatrix newb = parray[i];

					for (SimpleMatrix prevB : prevBs) {
						double num =
								prevB.transpose().mult(parray[i]).get(0) /
										prevB.transpose().mult(prevB).get(0);

						newb = newb.minus(prevB.scale(num));
					}

					barray[i] = newb;
				}
			}

			SimpleMatrix B =
					new SimpleMatrix(NOcc * NVirt, prevBs.size());
			SimpleMatrix P =
					new SimpleMatrix(NOcc * NVirt, prevBs.size());

			for (int i = 0; i < prevBs.size(); i++) {
				B.setColumn(i, 0, prevBs.get(i).getDDRM().data);

				P.setColumn(i, 0,
						prevBs.get(i).minus(prevPs.get(i)).getDDRM().data);
			}


			SimpleMatrix lhs = B.transpose().mult(P);

			SimpleMatrix rhs = B.transpose().mult(F);

			SimpleMatrix alpha = lhs.solve(rhs);

			for (int a = 0; a < xarray.length; a++) {
				rarray[a] = new SimpleMatrix(NOcc * NVirt, 1);
				xarray[a] = new SimpleMatrix(NOcc * NVirt, 1);
			}

			for (int i = 0; i < alpha.numRows(); i++) {
				for (int j = 0; j < alpha.numCols(); j++) {
					// B with tilde
					rarray[j] = rarray[j].plus((prevBs.get(i)
							.minus(prevPs.get(i))).scale(
							alpha.get(i, j)));
					xarray[j] = xarray[j].plus(
							prevBs.get(i).scale(alpha.get(i, j)));
				}
			}

			for (int j = 0; j < alpha.numCols(); j++) {

				// B0 is Farray, no tilde
				rarray[j] = rarray[j].minus(Farray[j]);

				xarray[j] = Dinv.mult(xarray[j]);

				if (mag(rarray[j]) < 1E-7) {
					iterable[j] = 1;
				}
				else if (Double.isNaN(mag(rarray[j]))) {
					System.err.println(
							"Pople algorithm fails; reverting to Thiel " +
									"algorithm (don't panic)...");
					return null;
				}
				else {
					iterable[j] = 0;
//					System.out.println("convergence test: " + mag(rarray[j]));
				}

			}

		}
		return xarray;
	}

	private static void testMain() throws InterruptedException {
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
		System.out.println(s);

		NanoStopWatch nsw = NanoStopWatch.sw();
		double time = 0;
		for (int i = 0; i < 1000; i++) {
			s = new SolutionR(atoms, rm).compute();
			nsw.start();
			s.compute();
			time += nsw.stop();
			TimeUnit.MILLISECONDS.sleep(1);
		}
//
		System.out.println("time = " + time / 1000);
//		System.out.println("s.energy = " + s.energy);
//		GeometryOptimization.of(s).compute();
//		System.out.println("s.energy = " + s.energy);
	}
}
