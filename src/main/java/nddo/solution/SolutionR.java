package nddo.solution;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.structs.MoleculeInfo;
import org.apache.logging.log4j.Level;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.Arrays;

import static nddo.State.config;
import static nddo.State.nom;
import static org.ejml.dense.row.CommonOps_DDRM.*;

public class SolutionR extends Solution {
	private final SimpleMatrix[] Farray, Darray, commutatorarray;
	private final double[] earray;
	private final SimpleMatrix B, Bforediis;
	private final transient DMatrixRMaj ddrm;

	public double[] integralArray;
	public SimpleMatrix F, E, Emat, C, COcc, CVirt, Ct, CtOcc, CtVirt;

	public SolutionR(MoleculeInfo mi, NDDOAtom[] atoms) {
		super(mi, atoms);

		Farray = new SimpleMatrix[LENGTH];
		Darray = new SimpleMatrix[LENGTH];
		commutatorarray = new SimpleMatrix[LENGTH];
		earray = new double[LENGTH];
		B = new SimpleMatrix(LENGTH, LENGTH);
		Bforediis = new SimpleMatrix(LENGTH, LENGTH);
		ddrm = new DMatrixRMaj(nOrbitals, nOrbitals);
	}

	@Override
	public Solution withNewAtoms(NDDOAtom[] newAtoms) {
		SolutionR s = new SolutionR(rm, newAtoms);

		s.E = E;
		s.Ct = Ct;

		s.densityMatrix = densityMatrix;

		s.compute();

		return s;
	}

	@Override
	public void computePrivate() {
		integralArray = new double[rm.nIntegrals];

		int intCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						integralArray[intCount] =
								nom.OneCenterERI(orbitals[j], orbitals[j], orbitals[l], orbitals[l]) -
										0.5 * nom.OneCenterERI(orbitals[j], orbitals[l], orbitals[j],
												orbitals[l]);
						intCount++;
					}

					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								integralArray[intCount] = nom.G(orbitals[j], orbitals[j], orbitals[l], orbitals[m]);
								intCount++;
							}
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					integralArray[intCount] =
							1.5 * nom.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]) -
									0.5 * nom.OneCenterERI(orbitals[j], orbitals[j], orbitals[k], orbitals[k]);
					intCount++;
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								integralArray[intCount] = nom.G(orbitals[j], orbitals[k], orbitals[l], orbitals[m]);
								intCount++;
							}
						}
					}
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							integralArray[intCount] = -0.5 * nom.G(orbitals[j], orbitals[l], orbitals[k], orbitals[m]);
							intCount++;
						}
					}
				}
			}
		}

		if (densityMatrix == null) {
			SimpleMatrix[] matrices = Utils.symEigen(H);

			E = matrices[1].diag();
			Ct = matrices[0].transposei();

			densityMatrix = calculateDensityMatrix();
		}

		SimpleMatrix G = new SimpleMatrix(nOrbitals, nOrbitals);

		double DIISError, threshold;
		int numIt = 0;

		while (true) {
			intCount = 0;

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							val += densityMatrix.get(l, l) * integralArray[intCount];
							intCount++;
						}

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									val += densityMatrix.get(l, m) * integralArray[intCount];
									intCount++;
								}
							}
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						val += densityMatrix.get(j, k) * integralArray[intCount];
						intCount++;

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									val += densityMatrix.get(l, m) * integralArray[intCount];
									intCount++;
								}
							}
						}
					}
					else {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
								val += densityMatrix.get(l, m) * integralArray[intCount];
								intCount++;
							}
						}
					}

					G.set(j, k, val);
					G.set(k, j, val);
				}
			}

			F = H.plus(G);


			if (numIt < LENGTH) {
				DIISError = updateDiis(numIt);
			}
			else {
				System.arraycopy(Farray, 1, Farray, 0, LENGTH1);
				System.arraycopy(Darray, 1, Darray, 0, LENGTH1);
				System.arraycopy(earray, 1, earray, 0, LENGTH1);
				System.arraycopy(commutatorarray, 1, commutatorarray, 0, LENGTH1);

				extract(B.getDDRM(), 1, LENGTH, 1, LENGTH, B.getDDRM(), 0, 0);
				B.setRow(LENGTH1, 0, ZEROS);
				B.setColumn(LENGTH1, 0, ZEROS);

				extract(Bforediis.getDDRM(), 1, LENGTH, 1, LENGTH, Bforediis.getDDRM(), 0, 0);
				Bforediis.setRow(LENGTH1, 0, ZEROS);
				Bforediis.setColumn(LENGTH1, 0, ZEROS);

				DIISError = updateDiis(LENGTH1);
			}


			int len = Math.min(LENGTH, numIt + 1);
			int len1 = len - 1;
			int ediisSize = len + 1;
			// if true do EDIIS else DIIS
			if (commutatorarray[len1].elementMax() > 0.01) {
				SimpleMatrix mat = new SimpleMatrix(ediisSize, ediisSize);

				for (int i = 0; i < len; i++) {
					for (int j = i; j < len; j++) {
						double v = Bforediis.get(i, j);
						mat.set(i, j, v);
						mat.set(j, i, v);
					}
				}

				double[] row = new double[ediisSize];
				double[] col = new double[ediisSize];
				Arrays.fill(row, 1);
				Arrays.fill(col, 1);
				mat.setColumn(len, 0, row);
				mat.setRow(len, 0, col);
				mat.set(len, len, 0);

				SimpleMatrix rhs = SimpleMatrix.ones(ediisSize, 1);
				for (int i = 0; i < len; i++) {
					rhs.set(i, earray[i]);
				}

				double bestE = 0;
				SimpleMatrix bestDIIS = null;
				for (int i = 0; i <= len1; i++) {
					for (int[] tbr : TBRS[i]) {
						try {
							SimpleMatrix newmat = removeElements(mat, tbr);
							SimpleMatrix newrhs = removeElements(rhs, tbr);
							SimpleMatrix tempEdiis = addRows(newmat.solve(newrhs), tbr);
							tempEdiis.set(len, 0);
							boolean nonNegative = !(elementMin(tempEdiis.getDDRM()) < 0);

							if (nonNegative) {
								double e = 0;

								for (int a = 0; a < len; a++) {
									e -= earray[a] * tempEdiis.get(a);

									for (int b = 0; b < len; b++) {
										e += 0.5 * tempEdiis.get(a) * tempEdiis.get(b) * Bforediis.get(a, b);
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

				final SimpleMatrix finalDIIS = bestDIIS;

				SimpleMatrix F = new SimpleMatrix(nOrbitals, nOrbitals);
				for (int i = 0; i < len; i++) {
					F.plusi(finalDIIS.get(i), Farray[i]);
				}

				SimpleMatrix[] matrices = Utils.symEigen(F);
				SimpleMatrix E = matrices[1].diag();
				SimpleMatrix Ct = matrices[0].transposei();

				if (!Ct.hasUncountable()) {
					this.E = E;
					this.Ct = Ct;
				}

				densityMatrix = calculateDensityMatrix();
			}
			else {
				SimpleMatrix mat = new SimpleMatrix(ediisSize, ediisSize);

				for (int i = 0; i < len; i++) {
					for (int j = i; j < len; j++) {
						mat.set(i, j, B.get(i, j));
						mat.set(j, i, B.get(i, j));

					}
				}

				double[] a = new double[ediisSize];
				Arrays.fill(a, 1);
				mat.setColumn(len, 0, a);
				mat.setRow(len, 0, a);
				mat.set(len, len, 0);

				SimpleMatrix rhs = new SimpleMatrix(ediisSize, 1);
				rhs.set(len, 0, 1);

				SimpleMatrix DIIS = mat.solve(rhs);

				SimpleMatrix F = new SimpleMatrix(nOrbitals, nOrbitals);
				SimpleMatrix D = new SimpleMatrix(nOrbitals, nOrbitals);

				for (int i = 0; i < len; i++) {
					double diis = DIIS.get(i);
					F.plusi(diis, Farray[i]);
					D.plusi(diis, Darray[i]);
				}

				SimpleMatrix[] matrices = Utils.symEigen(F);
				SimpleMatrix E = matrices[1].diag();
				SimpleMatrix Ct = matrices[0].transposei();

				if (!Ct.hasUncountable()) {
					this.E = E;
					this.Ct = Ct;
				}

				densityMatrix = calculateDensityMatrix();
			}

			if (numIt > config.rhf_numIt_max) {
				IllegalStateException e = new IllegalStateException(this.rm.debugName() + " unstable");
				rm.getLogger().error(e);
				throw e;
			}

			threshold = THRESHOLDS[numIt];
			boolean b = DIISError < threshold;

			Level level = b ? Level.DEBUG : Level.TRACE;
			rm.getLogger().log(level, "numIt: {}, DIISError: {}, threshold: {}", numIt, DIISError, threshold);

			if (b) break;

			numIt++;
		}
	}

	protected void findMatrices() {
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
	}

	protected void findEnergyAndHf() {
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

	protected void findHomoLumo() {
		homo = nElectrons > 0 ? E.get(rm.nOccAlpha - 1, 0) : 0;
		lumo = nElectrons != nOrbitals << 1 ? E.get(rm.nOccAlpha, 0) : 0;
	}

	private double updateDiis(int len1) {
		Farray[len1] = F;
		Darray[len1] = densityMatrix;
		earray[len1] = -0.5 * H.mult(densityMatrix).trace();

		commutatorarray[len1] = commutator(F, densityMatrix);
		double DIISError = commutatorarray[len1].normF();

		for (int i = 0; i <= len1; i++) {
			multTransB(commutatorarray[len1].getDDRM(), commutatorarray[i].getDDRM(), ddrm);
			double product = trace(ddrm);

			B.set(i, len1, product);
			B.set(len1, i, product);

			mult(Farray[i].getDDRM(), Darray[len1].getDDRM(), ddrm);
			double t1 = trace(ddrm);
			mult(Farray[len1].getDDRM(), Darray[i].getDDRM(), ddrm);
			double t2 = trace(ddrm);

			product = 0.5 * (t1 + t2);

			Bforediis.set(i, len1, product);
			Bforediis.set(len1, i, product);
		}

		return DIISError;
	}

	private SimpleMatrix calculateDensityMatrix() {
		CtOcc = Ct.extractMatrix(0, rm.nOccAlpha, 0, Ct.numCols());
		SimpleMatrix output = new SimpleMatrix(nOrbitals, nOrbitals);
		multInner(CtOcc.getDDRM(), output.getDDRM());
		return output.scalei(2);
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
