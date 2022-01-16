package nddo.solution;

import nddo.Constants;
import nddo.NDDOAtom;
import nddo.structs.MoleculeInfo;
import org.apache.logging.log4j.Level;
import org.ejml.data.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import static nddo.State.config;
import static nddo.State.nom;
import static org.ejml.dense.row.CommonOps_DDRM.*;

public class SolutionU extends Solution {
	private final SimpleMatrix[] Farrayalpha, Farraybeta, Darrayalpha, Darraybeta, commutatorarrayalpha,
			commutatorarraybeta;

	public SimpleMatrix Ea, Eb, Eamat, Ebmat, Fa, Fb,
			Ca, CaOcc, CaVirt, Cta, CtaOcc, CtaVirt, Cb, CbOcc, CbVirt, Ctb, CtbOcc, CtbVirt;
	public double[] integralArrayCoulomb, integralArrayExchange;

	protected SolutionU(MoleculeInfo mi, NDDOAtom[] atoms) {
		super(mi, atoms);

		if (nElectrons % 2 == mult % 2 || mult < 1) {
			mi.getLogger().error("Please check mult and charge: nElectrons: {}, mult: {}", nElectrons, mult);
		}

		Farrayalpha = new SimpleMatrix[8];
		Farraybeta = new SimpleMatrix[8];
		Darrayalpha = new SimpleMatrix[8];
		Darraybeta = new SimpleMatrix[8];
		commutatorarrayalpha = new SimpleMatrix[8];
		commutatorarraybeta = new SimpleMatrix[8];
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

		int intCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						integralArrayCoulomb[intCount] =
								nom.OneCenterERI(orbitals[j], orbitals[j], orbitals[l], orbitals[l]);
						intCount++;
					}
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								integralArrayCoulomb[intCount] =
										nom.G(orbitals[j], orbitals[j], orbitals[l], orbitals[m]);
								intCount++;
							}
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					integralArrayCoulomb[intCount] =
							2 * nom.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]);
					intCount++;
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						for (int m : missingOfAtom[atomOfOrb[j]]) {
							if (atomOfOrb[l] == atomOfOrb[m]) {
								integralArrayCoulomb[intCount] =
										nom.G(orbitals[j], orbitals[k], orbitals[l], orbitals[m]);
								intCount++;
							}
						}
					}
				}
			}
		}

		integralArrayExchange = new double[rm.nExchangeInts];

		intCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						integralArrayExchange[intCount] =
								-nom.OneCenterERI(orbitals[j], orbitals[l], orbitals[j], orbitals[l]);
						intCount++;
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					integralArrayExchange[intCount] =
							-nom.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]) -
									nom.OneCenterERI(orbitals[j], orbitals[j], orbitals[k], orbitals[k]);
					intCount++;
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						for (int m : orbsOfAtom[atomOfOrb[k]]) {
							integralArrayExchange[intCount] =
									-nom.G(orbitals[j], orbitals[l], orbitals[k], orbitals[m]);
							intCount++;
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

			alphaDensity = calculateDensityMatrix(Cta, rm.nOccAlpha);
			betaDensity = calculateDensityMatrix(Ctb, rm.nOccBeta);
		}

		SimpleMatrix J = new SimpleMatrix(nOrbitals, nOrbitals);
		SimpleMatrix Ka = new SimpleMatrix(nOrbitals, nOrbitals);
		SimpleMatrix Kb = new SimpleMatrix(nOrbitals, nOrbitals);

		double DIISError, threshold;
		int numIt = 0;

		while (true) {
			intCount = 0;

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							val += (alphaDensity.get(l, l) + betaDensity.get(l, l)) * integralArrayCoulomb[intCount];
							intCount++;
						}

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									val += (alphaDensity.get(l, m) + betaDensity.get(l, m)) *
											integralArrayCoulomb[intCount];
									intCount++;
								}
							}
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						val += (alphaDensity.get(j, k) + betaDensity.get(j, k)) * integralArrayCoulomb[intCount];
						intCount++;

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									val += (alphaDensity.get(l, m) + betaDensity.get(l, m)) *
											integralArrayCoulomb[intCount];
									intCount++;
								}
							}
						}
					}

					J.set(j, k, val);
					J.set(k, j, val);
				}
			}

			intCount = 0;

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double vala = 0;
					double valb = 0;
					if (j == k) {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							vala += alphaDensity.get(l, l) * integralArrayExchange[intCount];
							valb += betaDensity.get(l, l) * integralArrayExchange[intCount];
							intCount++;
						}

					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						vala += alphaDensity.get(j, k) * integralArrayExchange[intCount];
						valb += betaDensity.get(j, k) * integralArrayExchange[intCount];
						intCount++;

					}
					else {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
								vala += alphaDensity.get(l, m) * integralArrayExchange[intCount];
								valb += betaDensity.get(l, m) * integralArrayExchange[intCount];
								intCount++;
							}
						}
					}

					Ka.set(j, k, vala);
					Ka.set(k, j, vala);
					Kb.set(j, k, valb);
					Kb.set(k, j, valb);
				}
			}

			Fa = H.plus(J).plusi(Ka);
			Fb = H.plus(J).plusi(Kb);


			if (numIt < LENGTH) {
				DIISError = updateDiis(numIt);
			}
			else {
				System.arraycopy(Farrayalpha, 1, Farrayalpha, 0, LENGTH1);
				System.arraycopy(Darrayalpha, 1, Darrayalpha, 0, LENGTH1);
				System.arraycopy(commutatorarrayalpha, 1, commutatorarrayalpha, 0, LENGTH1);

				System.arraycopy(Farraybeta, 1, Farraybeta, 0, LENGTH1);
				System.arraycopy(Darraybeta, 1, Darraybeta, 0, LENGTH1);
				System.arraycopy(commutatorarraybeta, 1, commutatorarraybeta, 0, LENGTH1);

				System.arraycopy(earray, 1, earray, 0, LENGTH1);

				extract(B.getDDRM(), 1, LENGTH, 1, LENGTH, B.getDDRM(), 0, 0);
				B.setRow(LENGTH1, 0, ZEROS);
				B.setColumn(LENGTH1, 0, ZEROS);

				extract(Bforediis.getDDRM(), 1, LENGTH, 1, LENGTH, Bforediis.getDDRM(), 0, 0);
				Bforediis.setRow(LENGTH1, 0, ZEROS);
				Bforediis.setColumn(LENGTH1, 0, ZEROS);

				DIISError = updateDiis(LENGTH1);
			}


			int len = Math.min(LENGTH, numIt + 1);
			if (commutatorarrayalpha[len - 1].elementMax() > ediisThreshold
					|| commutatorarraybeta[len - 1].elementMax() > ediisThreshold) {
				fromDiis(ediis(len));
			}
			else {
				try {
					fromDiis(diis(len));
				} catch (SingularMatrixException e) {
					rm.getLogger().warn("DIIS has failed; using damping instead...");

					SimpleMatrix[] matricesa = Utils.symEigen(Fa);
					SimpleMatrix[] matricesb = Utils.symEigen(Fb);

					Ea = matricesa[1].diag();
					Eb = matricesb[1].diag();

					Cta = matricesa[0].transposei();
					Ctb = matricesb[0].transposei();

					double damp = config.uhf_damp;
					alphaDensity = calculateDensityMatrix(Cta, rm.nOccAlpha).scalei(1 - damp)
							.plusi(damp, alphaDensity);
					betaDensity = calculateDensityMatrix(Ctb, rm.nOccBeta).scalei(1 - damp)
							.plusi(damp, betaDensity);
				}
			}


			boolean exceededLimitButOK = false;
			if (numIt > config.uhf_numIt_max) {
				if (DIISError < config.diiserror_hard_limit) {
					rm.getLogger().warn("numIt limit reached but DIISError ({}) below hard limit;" +
									" finished prematurely.", DIISError);
					exceededLimitButOK = true;
				}
				else {
					IllegalStateException e = new IllegalStateException(this.rm.debugName() + " unstable");
					rm.getLogger().error("DIISError: {}, numIt: {}", DIISError, numIt, e);
					throw e;
				}
			}

			threshold = THRESHOLDS[numIt];
			boolean b = exceededLimitButOK || DIISError < threshold;

			Level level = b ? Level.DEBUG : Level.TRACE;
			rm.getLogger().log(level, "numIt: {}, DIISError: {}, threshold: {}", numIt, DIISError, threshold);

			if (b) break;

			numIt++;
		}
	}

	private SimpleMatrix calculateDensityMatrix(SimpleMatrix Ct, int nElectrons) {
		SimpleMatrix CtOcc = Ct.extractMatrix(0, nElectrons, 0, Ct.numCols());
		SimpleMatrix output = new SimpleMatrix(nOrbitals, nOrbitals);

		if (CtOcc.getNumElements() == 0) return output;

		multInner(CtOcc.getDDRM(), output.getDDRM());
		return output;
	}

	private double updateDiis(int len1) {
		Farrayalpha[len1] = Fa;
		Farraybeta[len1] = Fb;
		Darrayalpha[len1] = alphaDensity;
		Darraybeta[len1] = betaDensity;
		earray[len1] = -0.5 * (H.mult(alphaDensity).trace() + H.mult(betaDensity).trace());

		commutatorarrayalpha[len1] = commutator(Fa, alphaDensity);
		commutatorarraybeta[len1] = commutator(Fb, betaDensity);
		double DIISError = commutatorarrayalpha[len1].normF() + commutatorarraybeta[len1].normF();

		for (int i = 0; i <= len1; i++) {
			double product = 0;
			multTransB(commutatorarrayalpha[len1].getDDRM(), commutatorarrayalpha[i].getDDRM(), ddrm);
			product += trace(ddrm);
			multTransB(commutatorarraybeta[len1].getDDRM(), commutatorarraybeta[i].getDDRM(), ddrm);
			product += trace(ddrm);

			B.set(i, len1, product);
			B.set(len1, i, product);

			mult(Farrayalpha[i].getDDRM(), Darrayalpha[len1].getDDRM(), ddrm);
			double t1 = trace(ddrm);
			mult(Farrayalpha[len1].getDDRM(), Darrayalpha[i].getDDRM(), ddrm);
			double t2 = trace(ddrm);
			mult(Farraybeta[i].getDDRM(), Darraybeta[len1].getDDRM(), ddrm);
			double t3 = trace(ddrm);
			mult(Farraybeta[len1].getDDRM(), Darraybeta[i].getDDRM(), ddrm);
			double t4 = trace(ddrm);

			product = 0.5 * (t1 + t2 + t3 + t4);

			Bforediis.set(i, len1, product);
			Bforediis.set(len1, i, product);
		}

		return DIISError;
	}

	private void fromDiis(SimpleMatrix DIIS) {
		SimpleMatrix Fa = new SimpleMatrix(nOrbitals, nOrbitals);
		SimpleMatrix Fb = new SimpleMatrix(nOrbitals, nOrbitals);

		for (int i = 0, len = DIIS.getNumElements() - 1; i < len; i++) {
			double diis = DIIS.get(i);

			Fa.plusi(diis, Farrayalpha[i]);
			Fb.plusi(diis, Farraybeta[i]);
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

		alphaDensity = calculateDensityMatrix(Cta, rm.nOccAlpha);
		betaDensity = calculateDensityMatrix(Ctb, rm.nOccBeta);
	}

	@Override
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

	@Override
	protected void findHf() {
		energy = 0;

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

	@Override
	protected void findHomo() {
		homo = nElectrons > 0 ? Ea.get(rm.nOccAlpha - 1, 0) : 0;
		lumo = nElectrons != nOrbitals << 1 ? Eb.get(rm.nOccBeta, 0) : 0;
	}

	@Override
	public SimpleMatrix densityMatrix() {
		if (densityMatrix == null) densityMatrix = this.alphaDensity.plus(this.betaDensity);
		return densityMatrix;
	}
}
