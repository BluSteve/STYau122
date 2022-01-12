package nddo.structs;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MoleculeInfo { // low level molecule info representation
	public final int index; // must be unique, preferably small
	public final String name;
	public final boolean restricted;
	public final boolean useEdiis;
	public final int charge, mult;
	public final int[] atomicNumbers;
	public final int nElectrons, nOccAlpha, nOccBeta, nVirtAlpha, nVirtBeta, nonvAlpha, nonvBeta, nOrbitals,
			nIntegrals, nCoulombInts, nExchangeInts;
	public final int[][] orbsOfAtom, missingOfAtom;
	public final int[] atomOfOrb;
	public final int[] mats; // molecule atom types
	public final int[][] mnps; // molecule needed params
	public double[][] densityMatrices, densityMatricesExp;

	private transient String debugName;
	private transient Logger logger;

	public MoleculeInfo(int index, String name, boolean restricted, boolean useEdiis, int charge, int mult,
						int[] atomicNumbers, int nElectrons, int nOccAlpha, int nOccBeta, int nVirtAlpha,
						int nVirtBeta, int nonvAlpha, int nonvBeta, int nOrbitals, int nIntegrals, int nCoulombInts,
						int nExchangeInts, int[][] orbsOfAtom, int[][] missingOfAtom, int[] atomOfOrb, int[] mats,
						int[][] mnps, double[][] densityMatrices, double[][] densityMatricesExp) {
		this.index = index;
		this.name = name;
		this.restricted = restricted;
		this.useEdiis = useEdiis;
		this.charge = charge;
		this.mult = mult;
		this.atomicNumbers = atomicNumbers;
		this.nElectrons = nElectrons;
		this.nOccAlpha = nOccAlpha;
		this.nOccBeta = nOccBeta;
		this.nVirtAlpha = nVirtAlpha;
		this.nVirtBeta = nVirtBeta;
		this.nonvAlpha = nonvAlpha;
		this.nonvBeta = nonvBeta;
		this.nOrbitals = nOrbitals;
		this.nIntegrals = nIntegrals;
		this.nCoulombInts = nCoulombInts;
		this.nExchangeInts = nExchangeInts;
		this.orbsOfAtom = orbsOfAtom;
		this.missingOfAtom = missingOfAtom;
		this.atomOfOrb = atomOfOrb;
		this.mats = mats;
		this.mnps = mnps;
		this.densityMatrices = densityMatrices;
		this.densityMatricesExp = densityMatricesExp;
	}

	protected MoleculeInfo(MoleculeInfo mi) {
		this.index = mi.index;
		this.name = mi.name;
		this.restricted = mi.restricted;
		this.useEdiis = mi.useEdiis;
		this.charge = mi.charge;
		this.mult = mi.mult;
		this.atomicNumbers = mi.atomicNumbers;
		this.nElectrons = mi.nElectrons;
		this.nOccAlpha = mi.nOccAlpha;
		this.nOccBeta = mi.nOccBeta;
		this.nVirtAlpha = mi.nVirtAlpha;
		this.nVirtBeta = mi.nVirtBeta;
		this.nonvAlpha = mi.nonvAlpha;
		this.nonvBeta = mi.nonvBeta;
		this.nOrbitals = mi.nOrbitals;
		this.nIntegrals = mi.nIntegrals;
		this.nCoulombInts = mi.nCoulombInts;
		this.nExchangeInts = mi.nExchangeInts;
		this.orbsOfAtom = mi.orbsOfAtom;
		this.missingOfAtom = mi.missingOfAtom;
		this.atomOfOrb = mi.atomOfOrb;
		this.mats = mi.mats;
		this.mnps = mi.mnps;
		this.densityMatrices = mi.densityMatrices;
		this.densityMatricesExp = mi.densityMatricesExp;
	}

	public String debugName() {
		if (debugName == null)
			debugName = String.format("%03d-%s_%d_%s", index, name, charge, restricted ? "RHF" : "UHF");

		return debugName;
	}

	public Logger getLogger() {
		if (logger == null) {
			logger = LogManager.getLogger(String.format("%-26s", debugName()));
		}

		return logger;
	}

	public static class MIBuilder {
		public int index;
		public String name;
		public boolean restricted = true, useEdiis = true;
		public int charge, mult;
		public int[] atomicNumbers;
		public int nElectrons, nOrbitals;
		public int[][] orbsOfAtom, missingOfAtom;
		public int[] atomOfOrb;
		public int[] mats;
		public int[][] mnps;
		public double[][] densityMatrices, densityMatricesExp;

		public int findNIntegrals() {
			int size = 0;

			for (int j = 0; j < nOrbitals; j++) {
				for (int k = j; k < nOrbitals; k++) {
					if (j == k) {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							size++;
						}

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									size++;
								}
							}
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						size++;

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									size++;
								}
							}
						}
					}
					else {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
								size++;
							}
						}
					}
				}
			}

			return size;
		}

		@SuppressWarnings("DuplicatedCode")
		public int findNCoulombInts() {
			int size = 0;

			for (int j = 0; j < nOrbitals; j++) {
				for (int k = j; k < nOrbitals; k++) {
					if (j == k) {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							size++;
						}
						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									size++;
								}
							}
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						size++;
						for (int l : missingOfAtom[atomOfOrb[j]]) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (atomOfOrb[l] == atomOfOrb[m]) {
									size++;
								}
							}
						}
					}
				}
			}

			return size;
		}

		public int findNExchangeInts() {
			int size = 0;

			for (int j = 0; j < nOrbitals; j++) {
				for (int k = j; k < nOrbitals; k++) {
					if (j == k) {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							size++;
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						size++;
					}
					else {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
								size++;
							}
						}
					}
				}
			}

			return size;
		}

		public MoleculeInfo build() {
			if (isValid()) {
				int nIntegrals = 0, nCoulombInts = 0, nExchangeInts = 0;

				if (restricted) {
					nIntegrals = findNIntegrals();
				}
				else {
					nCoulombInts = findNCoulombInts();
					nExchangeInts = findNExchangeInts();
				}

				int nOccAlpha = (nElectrons - mult + 1) / 2 + mult - 1;
				int nOccBeta = nElectrons - nOccAlpha;
				int nVirtAlpha = nOrbitals - nOccAlpha;
				int nVirtBeta = nOrbitals - nOccBeta;
				int nonvAlpha = nOccAlpha * nVirtAlpha;
				int nonvBeta = nOccBeta * nVirtBeta;

				return new MoleculeInfo(index, name, restricted, useEdiis, charge, mult, atomicNumbers, nElectrons,
						nOccAlpha, nOccBeta, nVirtAlpha, nVirtBeta, nonvAlpha, nonvBeta, nOrbitals, nIntegrals,
						nCoulombInts, nExchangeInts, orbsOfAtom, missingOfAtom, atomOfOrb, mats, mnps,
						densityMatrices, densityMatricesExp);
			}

			throw new IllegalStateException("Invalid MIBuilder parameters for building!");
		}

		@SuppressWarnings("RedundantIfStatement")
		private boolean isValid() {
			if (mats.length != mnps.length) return false;

			int nAtoms = atomicNumbers.length;
			if (nAtoms != orbsOfAtom.length || nAtoms != missingOfAtom.length) return false;

			if (nOrbitals != atomOfOrb.length) return false;

			if (!restricted) {
				if (nElectrons % 2 == mult % 2 || mult < 1) {
					return false;
				}
			}

			return true;
		}
	}
}
