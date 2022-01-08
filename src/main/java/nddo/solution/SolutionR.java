package nddo.solution;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.structs.MoleculeInfo;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.Arrays;

import static nddo.State.nom;
import static nddo.State.config;

public class SolutionR extends Solution {

	public double[] integralArray;
	public SimpleMatrix C, COcc, CVirt, Ct, CtOcc, CtVirt, F, E, Emat;

	public SolutionR(MoleculeInfo mi, NDDOAtom[] atoms) {
		super(mi, atoms);
	}

	@Override
	public Solution withNewAtoms(NDDOAtom[] newAtoms) {
		return new SolutionR(rm, newAtoms).compute();
	}

	@Override
	public SolutionR compute() {
		integralArray = new double[rm.nIntegrals];

		int integralcount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) { // case 1
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						integralArray[integralcount] =
								nom.OneCenterERI(orbitals[j], orbitals[j], orbitals[l], orbitals[l]) -
										0.5 * nom.OneCenterERI(orbitals[j], orbitals[l], orbitals[j],
												orbitals[l]);
						integralcount++;
					}

					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								integralArray[integralcount] =
										nom.G(orbitals[j], orbitals[j], orbitals[l], orbitals[m]);
								integralcount++;
							}
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) { // case 2
					integralArray[integralcount] =
							1.5 * nom.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]) -
									0.5 * nom.OneCenterERI(orbitals[j], orbitals[j], orbitals[k], orbitals[k]);
					integralcount++;
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								integralArray[integralcount] =
										nom.G(orbitals[j], orbitals[k], orbitals[l], orbitals[m]);
								integralcount++;
							}
						}
					}
				}
				else { // case 3
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							integralArray[integralcount] = -0.5 * nom.G(orbitals[j], orbitals[l],
									orbitals[k], orbitals[m]);
							integralcount++;
						}
					}
				}
			}
		}

		SimpleMatrix[] matrices = Utils.symEigen(H);
		E = matrices[1].diag();
		Ct = matrices[0].transpose();
		SimpleMatrix G = new SimpleMatrix(Ct.numRows(), Ct.numCols());
		F = H.copy();

		densityMatrix = calculateDensityMatrix(Ct);
		SimpleMatrix olddensity;

		SimpleMatrix[] Farray = new SimpleMatrix[8];
		SimpleMatrix[] Darray = new SimpleMatrix[8];

		double[] earray = new double[8];

		SimpleMatrix Bforediis = new SimpleMatrix(8, 8);
		SimpleMatrix B = new SimpleMatrix(8, 8);

		SimpleMatrix[] commutatorarray = new SimpleMatrix[8];

		int numIt = 0;
		double DIISError = 10;
		int itSinceLastDIIS = 0;

		while (true) {
			double threshold = sigmoid(numIt);
			if (!(DIISError > threshold)) break;

			olddensity = densityMatrix;
			integralcount = 0;

			// this entire block of code fills up the G matrix, and it calls
			// the integralarray to save time.
			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							val += densityMatrix.get(l, l) * integralArray[integralcount];
							integralcount++;
						}

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									val += densityMatrix.get(l, m) * integralArray[integralcount];
									integralcount++;
								}
							}
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						val += densityMatrix.get(j, k) * integralArray[integralcount];
						integralcount++;

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									val += densityMatrix.get(l, m) * integralArray[integralcount];
									integralcount++;
								}
							}
						}
					}
					else {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
								val += densityMatrix.get(l, m) * integralArray[integralcount];
								integralcount++;
							}
						}
					}

					G.set(j, k, val);
					G.set(k, j, val);
				}
			}

			F = H.plus(G);

			if (numIt < Farray.length) {
				Farray[numIt] = F;

				Darray[numIt] = densityMatrix;
				earray[numIt] = -0.5 * H.mult(densityMatrix).diag().elementSum();
				commutatorarray[numIt] = commutator(F, densityMatrix);
				DIISError = commutatorarray[numIt].normF();

				for (int i = 0; i <= numIt; i++) {
					double product = commutatorarray[numIt].mult(commutatorarray[i].transpose()).diag().elementSum();

					B.set(i, numIt, product);
					B.set(numIt, i, product);

					product = 0.5 * Farray[i].mult(Darray[numIt]).diag().elementSum() +
							Farray[numIt].mult(Darray[i]).diag().elementSum();

					Bforediis.set(i, numIt, product);
					Bforediis.set(numIt, i, product);
				}
			}
			else {
				for (int i = 0; i < Farray.length - 1; i++) {
					Farray[i] = Farray[i + 1];
					Darray[i] = Darray[i + 1];
					commutatorarray[i] = commutatorarray[i + 1];
					earray[i] = earray[i + 1];
				}

				Farray[Farray.length - 1] = F;
				Darray[Darray.length - 1] = densityMatrix;
				earray[Darray.length - 1] = -0.5 * H.mult(densityMatrix).diag().elementSum();

				commutatorarray[Darray.length - 1] = commutator(F, densityMatrix);
				DIISError = commutatorarray[Darray.length - 1].normF();

				// B is dy/dx sort of, make dy/dx 0
				SimpleMatrix newB = new SimpleMatrix(8, 8);
				SimpleMatrix newBforediis = new SimpleMatrix(8, 8);

				for (int i = 0; i < Farray.length - 1; i++) {
					for (int j = i; j < Farray.length - 1; j++) {
						newB.set(i, j, B.get(i + 1, j + 1));
						newB.set(j, i, B.get(i + 1, j + 1));
						newBforediis.set(i, j, Bforediis.get(i + 1, j + 1));
						newBforediis.set(j, i, Bforediis.get(i + 1, j + 1));
					}
				}

				for (int i = 0; i < Farray.length; i++) {
					double product = commutatorarray[Farray.length - 1].transpose().mult(commutatorarray[i])
							.diag().elementSum();

					newB.set(i, Farray.length - 1, product);
					newB.set(Farray.length - 1, i, product);

					product = 0.5 * Farray[i].mult(Darray[Farray.length - 1]).diag().elementSum() +
							Farray[Farray.length - 1].mult(Darray[i]).diag().elementSum();

					newBforediis.set(i, Farray.length - 1, product);
					newBforediis.set(Farray.length - 1, i, product);
				}

				B = newB;
				Bforediis = newBforediis;
			}

			int ediisSize = Math.min(Farray.length + 1, numIt + 2);

			// if true do EDIIS else DIIS
			if (commutatorarray[Math.min(Farray.length - 1, numIt)].elementMax() > 0.01 && itSinceLastDIIS < 50) {
				itSinceLastDIIS++;
				SimpleMatrix mat = new SimpleMatrix(ediisSize, ediisSize);

				for (int i = 0; i < ediisSize - 1; i++) {
					for (int j = i; j < ediisSize - 1; j++) {
						mat.set(i, j, Bforediis.get(i, j));
						mat.set(j, i, Bforediis.get(i, j));
					}
				}

				double[] row = new double[mat.numRows()];
				double[] col = new double[mat.numCols()];
				Arrays.fill(row, 1);
				Arrays.fill(col, 1);
				mat.setColumn(mat.numCols() - 1, 0, row);
				mat.setRow(mat.numRows() - 1, 0, col);
				mat.set(mat.numRows() - 1, mat.numCols() - 1, 0);

				SimpleMatrix rhs = SimpleMatrix.ones(mat.numRows(), 1);
				for (int i = 0; i < ediisSize - 1; i++) {
					rhs.set(i, earray[i]);
				}

				double bestE = 0;
				SimpleMatrix bestDIIS = null;
				int n = mat.numRows() - 2;
				for (int i = 0; i <= n; i++) {
					for (int[] tbr : TBRS[i]) {
						try {
							SimpleMatrix newmat = removeElements(mat, tbr);
							SimpleMatrix newrhs = removeElements(rhs, tbr);
							SimpleMatrix tempEdiis = addRows(newmat.solve(newrhs), tbr);
							tempEdiis.set(tempEdiis.numRows() - 1, 0);
							boolean nonNegative = !(CommonOps_DDRM.elementMin(tempEdiis.getDDRM()) < 0);

							if (nonNegative) {
								double e = 0;

								for (int a = 0; a < tempEdiis.getNumElements() - 1; a++) {
									e -= earray[a] * tempEdiis.get(a);
								}

								for (int a = 0; a < tempEdiis.getNumElements() - 1; a++) {
									for (int b = 0; b < tempEdiis.getNumElements() - 1; b++) {
										e += 0.5 * tempEdiis.get(a) * tempEdiis.get(b) * 0.5 * B.get(a, b);
									}
								}

								if (e < bestE) {
									bestE = e;
									bestDIIS = tempEdiis;
								}
							}

						} catch (SingularMatrixException ignored) {
						}
					}
				}

				SimpleMatrix finalDIIS = bestDIIS;

				SimpleMatrix F = new SimpleMatrix(densityMatrix.numRows(), densityMatrix.numCols());

				for (int i = 0; i < finalDIIS.getNumElements() - 1; i++) {
					F = F.plus(Farray[i].scale(finalDIIS.get(i)));
				}

				matrices = Utils.symEigen(F);
				E = matrices[1].diag();
				Ct = matrices[0].transpose();

				if (Ct.get(0, 0) != Ct.get(0, 0)) {
					matrices = Utils.symEigen(this.F);
					E = matrices[1].diag();
					Ct = matrices[0].transpose();
				}

				densityMatrix = calculateDensityMatrix(Ct);
			}
			else {
				itSinceLastDIIS = 0;
				SimpleMatrix mat = new SimpleMatrix(ediisSize, ediisSize);

				for (int i = 0; i < Math.min(Farray.length, numIt + 1); i++) {
					for (int j = i; j < Math.min(Farray.length, numIt + 1); j++) {
						mat.set(i, j, B.get(i, j));
						mat.set(j, i, B.get(i, j));

					}
				}

				double[] a = new double[mat.numRows()];
				Arrays.fill(a, 1);
				mat.setColumn(mat.numCols() - 1, 0, a);
				mat.setRow(mat.numRows() - 1, 0, a);
				mat.set(mat.numRows() - 1, mat.numCols() - 1, 0);

				SimpleMatrix rhs = new SimpleMatrix(mat.numRows(), 1);
				rhs.set(mat.numRows() - 1, 0, 1);

				try {
					SimpleMatrix DIIS = mat.solve(rhs);

					SimpleMatrix F = new SimpleMatrix(densityMatrix.numRows(), densityMatrix.numCols());

					SimpleMatrix D = new SimpleMatrix(densityMatrix.numRows(), densityMatrix.numCols());

					for (int i = 0; i < DIIS.getNumElements() - 1; i++) {
						F.plusi(Farray[i].scale(DIIS.get(i)));
						D.plusi(Darray[i].scale(DIIS.get(i)));
					}

					matrices = Utils.symEigen(F);
					E = matrices[1].diag();
					Ct = matrices[0].transpose();

					if (Ct.get(0, 0) != Ct.get(0, 0)) {
						matrices = Utils.symEigen(this.F);
						E = matrices[1].diag();
						Ct = matrices[0].transpose();
					}

					densityMatrix = calculateDensityMatrix(Ct);
				} catch (SingularMatrixException e) { // todo fix adrian
					matrices = Utils.symEigen(F);
					E = matrices[1].diag();
					Ct = matrices[0].transpose();

					double damp = 0.8;
					densityMatrix = calculateDensityMatrix(Ct).scale(1 - damp).plus(olddensity.scale(damp));
				}
			}

			if (numIt > config.rhf_numIt_max) {
				IllegalStateException e = new IllegalStateException(this.rm.debugName() + " unstable");
				rm.getLogger().error(e);
				throw e;
			}

			rm.getLogger().trace("numIt: {}, DIISError: {}, threshold: {}", numIt, DIISError, threshold);
			numIt++;
		}

		// todo make these precomputations optional
		CtOcc = Ct.extractMatrix(0, rm.nOccAlpha, 0, Ct.numCols());
		CtVirt = Ct.extractMatrix(rm.nOccAlpha, Ct.numCols(), 0, Ct.numCols());
		C = Ct.transpose();
		COcc = CtOcc.transpose();
		CVirt = CtVirt.transpose();

		SimpleMatrix sm = E.extractMatrix(rm.nOccAlpha, nOrbitals, 0, 1);
		sm.reshape(1, rm.nVirtAlpha);

		Emat = new SimpleMatrix(rm.nOccAlpha, rm.nVirtAlpha);
		for (int i = 0; i < rm.nOccAlpha; i++) {
			Emat.insertIntoThis(i, 0, sm.minus(E.get(i)));
		}

		findEnergyAndHf();
		if (nElectrons > 0) homo = E.get(nElectrons / 2 - 1, 0);
		else homo = 0;
		if (nElectrons != nOrbitals << 1) lumo = E.get(nElectrons / 2, 0);
		else lumo = 0;
		findDipole();

		return this;
	}

	private void findEnergyAndHf() {
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				energy += 0.5 * densityMatrix.get(j, k) * (H.get(j, k) + F.get(j, k));
			}
		}

		double heat = 0;
		for (int j = 0; j < atoms.length; j++) {
			heat += atoms[j].getAtomProperties().getHeat() - atoms[j].getParams().getEisol();
			for (int k = j + 1; k < atoms.length; k++) {
				energy += atoms[j].crf(atoms[k]);
			}
		}

		heat += energy;
		hf = heat / Constants.HEATCONV;
	}

	private SimpleMatrix calculateDensityMatrix(SimpleMatrix c) {
		SimpleMatrix densityMatrix = new SimpleMatrix(orbitals.length, orbitals.length);

		for (int i = 0; i < orbitals.length; i++) {
			for (int j = 0; j < orbitals.length; j++) {
				double sum = 0;
				int count = nElectrons;
				int counter = -1;
				while (count > 0) {
					counter++;
					sum += 2 * c.get(counter, i) * c.get(counter, j);
					count -= 2;
				}
				densityMatrix.set(i, j, sum);
			}
		}

		return densityMatrix;
	}

	@Override
	public SimpleMatrix alphaDensity() {
		if (alphaDensity == null) alphaDensity = densityMatrix.scale(0.5);
		return super.alphaDensity();
	}

	@Override
	public SimpleMatrix betaDensity() {
		if (betaDensity == null) betaDensity = densityMatrix.scale(0.5);
		return super.betaDensity();
	}
}
