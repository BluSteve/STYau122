package nddo.solution;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.structs.MoleculeInfo;
import org.apache.logging.log4j.Level;
import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.Arrays;

import static nddo.State.config;
import static nddo.State.nom;
import static org.ejml.dense.row.CommonOps_DDRM.multInner;


public class SolutionU extends Solution {
	private final int nAlpha = rm.nOccAlpha;
	private final int nBeta = rm.nOccBeta;
	public SimpleMatrix Ea, Eb, Eamat, Ebmat, Fa, Fb,
			Ca, CaOcc, CaVirt, Cta, CtaOcc, CtaVirt, Cb, CbOcc, CbVirt, Ctb, CtbOcc, CtbVirt;
	public double[] integralArrayCoulomb, integralArrayExchange;

	protected SolutionU(MoleculeInfo mi, NDDOAtom[] atoms) {
		super(mi, atoms);
		if (nElectrons % 2 == mult % 2 || mult < 1) {
			mi.getLogger().error("Please check mult and charge: nElectrons: {}, mult: {}", nElectrons, mult);
		}
	}

	@Override
	public Solution withNewAtoms(NDDOAtom[] newAtoms) {
		SolutionU s = new SolutionU(rm, newAtoms);

		s.Ea = Ea;
		s.Eb = Eb;
		s.Cta = Cta;
		s.Ctb = Ctb;

		s.alphaDensity = alphaDensity;
		s.betaDensity = betaDensity;

		s.compute();

		return s;
	}

	@Override
	public void computePrivate() {
		integralArrayCoulomb = new double[rm.nCoulombInts];

		int integralCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						integralArrayCoulomb[integralCount] =
								nom.OneCenterERI(orbitals[j], orbitals[j], orbitals[l], orbitals[l]);
						integralCount++;
					}
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								integralArrayCoulomb[integralCount] =
										nom.G(orbitals[j], orbitals[j], orbitals[l], orbitals[m]);
								integralCount++;
							}
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					integralArrayCoulomb[integralCount] =
							2 * nom.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]);
					integralCount++;
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								integralArrayCoulomb[integralCount] =
										nom.G(orbitals[j], orbitals[k], orbitals[l], orbitals[m]);
								integralCount++;
							}
						}
					}
				}
			}
		}

		integralArrayExchange = new double[rm.nExchangeInts];

		integralCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						integralArrayExchange[integralCount] =
								-nom.OneCenterERI(orbitals[j], orbitals[l], orbitals[j], orbitals[l]);
						integralCount++;
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					integralArrayExchange[integralCount] =
							-nom.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]) -
									nom.OneCenterERI(orbitals[j], orbitals[j], orbitals[k], orbitals[k]);
					integralCount++;
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							integralArrayExchange[integralCount] =
									-nom.G(orbitals[j], orbitals[l], orbitals[k], orbitals[m]);
							integralCount++;
						}
					}
				}
			}
		}

		if (alphaDensity == null || betaDensity == null) {
			SimpleMatrix[] matrices = Utils.symEigen(H);

			Ea = matrices[1].diag();
			Eb = Ea.copy();

			Cta = matrices[0].transposei();
			Ctb = Cta.copy();

			alphaDensity = calculateDensityMatrix(Cta, nAlpha);
			betaDensity = calculateDensityMatrix(Ctb, nBeta);
		}

		SimpleMatrix j1 = new SimpleMatrix(nOrbitals, nOrbitals);
		SimpleMatrix ka = new SimpleMatrix(nOrbitals, nOrbitals);
		SimpleMatrix kb = new SimpleMatrix(nOrbitals, nOrbitals);
		SimpleMatrix oldalphadensity, oldbetadensity;

		int Jcount, Kcount;

		int numIt = 0;

		SimpleMatrix[] Farrayalpha = new SimpleMatrix[8];
		SimpleMatrix[] Farraybeta = new SimpleMatrix[8];

		SimpleMatrix[] Darrayalpha = new SimpleMatrix[8];
		SimpleMatrix[] Darraybeta = new SimpleMatrix[8];


		SimpleMatrix B = new SimpleMatrix(8, 8);
		SimpleMatrix Bforediis = new SimpleMatrix(8, 8);


		SimpleMatrix[] commutatorarrayalpha = new SimpleMatrix[8];
		SimpleMatrix[] commutatorarraybeta = new SimpleMatrix[8];
		double[] earray = new double[8];

		double DIISError, threshold;

		while (true) {
			oldalphadensity = alphaDensity.copy();
			oldbetadensity = betaDensity.copy();

			Jcount = 0;
			Kcount = 0;

			//construct J matrix

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {

						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							val += (alphaDensity.get(l, l) +
									betaDensity.get(l, l)) *
									integralArrayCoulomb[Jcount];
							Jcount++;
						}

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									val += (alphaDensity.get(l, m) +
											betaDensity.get(l, m)) *
											integralArrayCoulomb[Jcount];
									Jcount++;
								}
							}
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						val += (alphaDensity.get(j, k) +
								betaDensity.get(j, k)) *
								integralArrayCoulomb[Jcount];
						Jcount++;

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									val += (alphaDensity.get(l, m) +
											betaDensity.get(l, m)) *
											integralArrayCoulomb[Jcount];
									Jcount++;
								}
							}
						}
					}


					j1.set(j, k, val);
					j1.set(k, j, val);
				}
			}

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double vala = 0;
					double valb = 0;
					if (j == k) {

						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							vala += alphaDensity.get(l, l) *
									integralArrayExchange[Kcount];
							valb += betaDensity.get(l, l) *
									integralArrayExchange[Kcount];
							Kcount++;
						}

					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						vala += alphaDensity.get(j, k) *
								integralArrayExchange[Kcount];
						valb += betaDensity.get(j, k) *
								integralArrayExchange[Kcount];
						Kcount++;

					}
					else {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
								vala += alphaDensity.get(l, m) *
										integralArrayExchange[Kcount];
								valb += betaDensity.get(l, m) *
										integralArrayExchange[Kcount];
								Kcount++;
							}
						}
					}

					ka.set(j, k, vala);
					ka.set(k, j, vala);
					kb.set(j, k, valb);
					kb.set(k, j, valb);
				}
			}

			Fa = H.plus(j1).plus(ka);
			Fb = H.plus(j1).plus(kb);

			if (numIt < Farrayalpha.length) {

				Farrayalpha[numIt] = Fa.copy();
				Farraybeta[numIt] = Fb.copy();


				Darrayalpha[numIt] = alphaDensity.copy();
				Darraybeta[numIt] = betaDensity.copy();
				earray[numIt] =
						-0.5 * (H.mult(alphaDensity).diag().elementSum() + H.mult(betaDensity).diag().elementSum());

				commutatorarrayalpha[numIt] = commutator(Fa.copy(), alphaDensity.copy());
				commutatorarraybeta[numIt] = commutator(Fb.copy(), betaDensity.copy());

				DIISError = commutatorarrayalpha[numIt].normF() + commutatorarraybeta[numIt].normF();

				for (int i = 0; i <= numIt; i++) {

					double product =
							commutatorarrayalpha[numIt].mult(commutatorarrayalpha[i].transpose()).diag().elementSum() +
									commutatorarraybeta[numIt].mult(commutatorarraybeta[i].transpose()).diag()
											.elementSum();
					B.set(i, numIt, product);
					B.set(numIt, i, product);

					product = 0.5 * (Farrayalpha[i].mult(Darrayalpha[numIt]).diag().elementSum() +
							Farraybeta[i].mult(Darraybeta[numIt]).diag().elementSum()
							+ Farrayalpha[numIt].mult(Darrayalpha[i]).diag().elementSum() +
							Farraybeta[numIt].mult(Darraybeta[i]).diag().elementSum());

					Bforediis.set(i, numIt, product);
					Bforediis.set(numIt, i, product);

				}
			}

			else {

				for (int i = 0; i < Farrayalpha.length - 1; i++) {

					Farrayalpha[i] = Farrayalpha[i + 1].copy();
					Farraybeta[i] = Farraybeta[i + 1].copy();

					Darrayalpha[i] = Darrayalpha[i + 1].copy();
					Darraybeta[i] = Darraybeta[i + 1].copy();

					commutatorarrayalpha[i] = commutatorarrayalpha[i + 1].copy();
					commutatorarraybeta[i] = commutatorarraybeta[i + 1].copy();

					earray[i] = earray[i + 1];
				}

				Farrayalpha[Farrayalpha.length - 1] = Fa.copy();
				Farraybeta[Farraybeta.length - 1] = Fb.copy();

				Darrayalpha[Darrayalpha.length - 1] = alphaDensity.copy();
				Darraybeta[Darraybeta.length - 1] = betaDensity.copy();

				commutatorarrayalpha[Darrayalpha.length - 1] = commutator(Fa.copy(), alphaDensity.copy());
				commutatorarraybeta[Darraybeta.length - 1] = commutator(Fb.copy(), betaDensity.copy());

				earray[Darrayalpha.length - 1] =
						-0.5 * (H.mult(alphaDensity).diag().elementSum() + H.mult(betaDensity).diag().elementSum());

				DIISError = commutatorarrayalpha[Darrayalpha.length - 1].normF() +
						commutatorarraybeta[Darraybeta.length - 1].normF();

				// B is dy/dx sort of, make dy/dx 0
				SimpleMatrix newB = new SimpleMatrix(8, 8);
				SimpleMatrix newBforediis = new SimpleMatrix(8, 8);

				for (int i = 0; i < Farrayalpha.length - 1; i++) {
					for (int j = i; j < Farrayalpha.length - 1; j++) {
						newB.set(i, j, B.get(i + 1, j + 1));
						newB.set(j, i, B.get(i + 1, j + 1));
						newBforediis.set(i, j, Bforediis.get(i + 1, j + 1));
						newBforediis.set(j, i, Bforediis.get(i + 1, j + 1));
					}
				}

				for (int i = 0; i < Farrayalpha.length; i++) {

					double product =
							commutatorarrayalpha[Farrayalpha.length - 1].transpose().mult(commutatorarrayalpha[i])
									.diag().elementSum() +
									commutatorarraybeta[Farraybeta.length - 1].transpose().mult(commutatorarraybeta[i])
											.diag().elementSum();

					newB.set(i, Farrayalpha.length - 1, product);
					newB.set(Farrayalpha.length - 1, i, product);

					product = 0.5 * (Farrayalpha[i].mult(Darrayalpha[Farrayalpha.length - 1]).diag().elementSum() +
							Farraybeta[i].mult(Darraybeta[Farrayalpha.length - 1]).diag().elementSum()
							+ Farrayalpha[Farrayalpha.length - 1].mult(Darrayalpha[i]).diag().elementSum() +
							Farraybeta[Farrayalpha.length - 1].mult(Darraybeta[i]).diag().elementSum());

					newBforediis.set(i, Farrayalpha.length - 1, product);
					newBforediis.set(Farrayalpha.length - 1, i, product);
				}

				B = newB.copy();
				Bforediis = newBforediis.copy();
			}

			int ediisSize = Math.min(Farrayalpha.length + 1, numIt + 2);

			if (commutatorarrayalpha[Math.min(Farrayalpha.length - 1, numIt)].elementMax() > 1E-3 ||
					commutatorarraybeta[Math.min(Farrayalpha.length - 1, numIt)].elementMax() > 1E-3) {


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


				SimpleMatrix finalDIIS = bestDIIS;

				SimpleMatrix Fa = new SimpleMatrix(alphaDensity.numRows(), alphaDensity.numCols());
				SimpleMatrix Fb = new SimpleMatrix(betaDensity.numRows(), betaDensity.numCols());


				SimpleMatrix Da = new SimpleMatrix(alphaDensity.numRows(), alphaDensity.numCols());
				SimpleMatrix Db = new SimpleMatrix(betaDensity.numRows(), betaDensity.numCols());

				for (int i = 0; i < finalDIIS.getNumElements() - 1; i++) {
					Fa.plusi(Farrayalpha[i].scale(finalDIIS.get(i)));
					Da.plusi(Darrayalpha[i].scale(finalDIIS.get(i)));

					Fb.plusi(Farraybeta[i].scale(finalDIIS.get(i)));
					Db.plusi(Darraybeta[i].scale(finalDIIS.get(i)));
				}

				SimpleMatrix[] matrices1 = Utils.symEigen(Fa);

				SimpleMatrix[] matrices2 = Utils.symEigen(Fb);

				Ea = matrices1[1].diag();

				Eb = matrices2[1].diag();

				Cta = matrices1[0].transpose();
				Ctb = matrices2[0].transpose();


				if (Cta.get(0, 0) != Cta.get(0, 0)) {

					matrices1 = Utils.symEigen(this.Fa);

					matrices2 = Utils.symEigen(Fb);

					Ea = matrices1[1].diag();

					Eb = matrices2[1].diag();

					Cta = matrices1[0].transpose();
					Ctb = matrices2[0].transpose();


				}

				alphaDensity = calculateDensityMatrix(Cta, nAlpha);

				betaDensity = calculateDensityMatrix(Ctb, nBeta);

			}

			else {
				SimpleMatrix mat = new SimpleMatrix(ediisSize, ediisSize);

				for (int i = 0; i < Math.min(Farrayalpha.length, numIt + 1); i++) {
					for (int j = i; j < Math.min(Farrayalpha.length, numIt + 1);
						 j++) {
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

				SimpleMatrix DIIS = mat.solve(rhs);

				SimpleMatrix Fa = new SimpleMatrix(nOrbitals, nOrbitals);
				SimpleMatrix Fb = new SimpleMatrix(nOrbitals, nOrbitals);

				SimpleMatrix Da = new SimpleMatrix(nOrbitals, nOrbitals);
				SimpleMatrix Db = new SimpleMatrix(nOrbitals, nOrbitals);

				for (int i = 0; i < DIIS.getNumElements() - 1; i++) {
					double diis = DIIS.get(i);
					Fa.plusi(diis, Farrayalpha[i]);
					Da.plusi(diis, Darrayalpha[i]);

					Fb.plusi(diis, Farraybeta[i]);
					Db.plusi(diis, Darraybeta[i]);
				}

				SimpleMatrix[] matrices1 = Utils.symEigen(Fa);
				SimpleMatrix[] matrices2 = Utils.symEigen(Fb);

				SimpleMatrix Ea = matrices1[1].diag();
				SimpleMatrix Eb = matrices2[1].diag();

				SimpleMatrix Cta = matrices1[0].transposei();
				SimpleMatrix Ctb = matrices2[0].transposei();

				if (!Cta.hasUncountable() && !Ctb.hasUncountable()) {
					this.Ea = Ea;
					this.Eb = Eb;
					this.Cta = Cta;
					this.Ctb = Ctb;
				}

				alphaDensity = calculateDensityMatrix(Cta, nAlpha);
				betaDensity = calculateDensityMatrix(Ctb, nBeta);
			}


			if (numIt > config.uhf_numIt_max) {
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
		CtaOcc = Cta.extractMatrix(0, rm.nOccAlpha, 0, Cta.numCols());
		CtaVirt = Cta.extractMatrix(rm.nOccAlpha, Cta.numCols(), 0, Cta.numCols());
		Ca = Cta.transpose();
		CaOcc = CtaOcc.transpose();
		CaVirt = CtaVirt.transpose();

		SimpleMatrix sm = Ea.extractMatrix(rm.nOccAlpha, nOrbitals, 0, 1);
		sm.reshape(1, rm.nVirtAlpha);

		Eamat = new SimpleMatrix(rm.nOccAlpha, rm.nVirtAlpha);
		for (int i = 0; i < rm.nOccAlpha; i++) {
			Eamat.insertIntoThis(i, 0, sm.minus(Ea.get(i)));
		}

		CtbOcc = Ctb.extractMatrix(0, rm.nOccBeta, 0, Ctb.numCols());
		CtbVirt = Ctb.extractMatrix(rm.nOccBeta, Ctb.numCols(), 0, Ctb.numCols());
		Cb = Ctb.transpose();
		CbOcc = CtbOcc.transpose();
		CbVirt = CtbVirt.transpose();

		sm = Eb.extractMatrix(rm.nOccBeta, nOrbitals, 0, 1);
		sm.reshape(1, rm.nVirtBeta);

		Ebmat = new SimpleMatrix(rm.nOccBeta, rm.nVirtBeta);
		for (int i = 0; i < rm.nOccBeta; i++) {
			Ebmat.insertIntoThis(i, 0, sm.minus(Eb.get(i)));
		}
	}

	protected void findEnergyAndHf() {
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				double v = H.get(j, k);
				energy += 0.5 * alphaDensity.get(j, k) * (v + Fa.get(j, k));
				energy += 0.5 * betaDensity.get(j, k) * (v + Fb.get(j, k));
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
		homo = nElectrons > 0 ? Ea.get(nAlpha - 1, 0) : 0;
		lumo = nElectrons != nOrbitals << 1 ? Eb.get(nBeta, 0) : 0;
	}

	private SimpleMatrix calculateDensityMatrix(SimpleMatrix Ct, int nElectrons) {
		SimpleMatrix CtOcc = Ct.extractMatrix(0, nElectrons, 0, Ct.numCols());
		SimpleMatrix output = new SimpleMatrix(nOrbitals, nOrbitals);
		multInner(CtOcc.getDDRM(), output.getDDRM());
		return output;
	}

	@Override
	public SimpleMatrix alphaDensity() {
		return this.alphaDensity;
	}

	@Override
	public SimpleMatrix betaDensity() {
		return this.betaDensity;
	}

	@Override
	public SimpleMatrix densityMatrix() {
		if (densityMatrix == null) densityMatrix = this.alphaDensity.plus(this.betaDensity);
		return densityMatrix;
	}
}
