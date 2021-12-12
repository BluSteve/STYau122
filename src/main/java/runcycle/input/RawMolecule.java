package runcycle.input;

import nddo.MoleculeInfo;
import scf.AtomHandler;
import scf.AtomProperties;
import scf.OrbitalProperties;
import tools.Utils;

import java.io.IOException;
import java.util.*;

public class RawMolecule extends MoleculeInfo { // for storage purposes
	public RawAtom[] atoms, expGeom;

	private RawMolecule() {
	}

	public static class RMBuilder {
		private int index;
		private boolean restricted = true;
		private int charge, mult;
		private double[] datum;
		private RawAtom[] atoms, expGeom;

		public static int findNIntegrals(MoleculeInfo mi) {
			int size = 0;

			for (int j = 0; j < mi.nOrbitals; j++) {
				for (int k = j; k < mi.nOrbitals; k++) {
					if (j == k) {
						for (int l : mi.orbsOfAtom[mi.atomOfOrb[j]]) {
							if (l > -1) {
								size++;
							}
						}

						for (int l : mi.missingOfAtom[mi.atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : mi.missingOfAtom[mi.atomOfOrb[j]]) {
									if (m > -1) {
										if (mi.atomOfOrb[l] == mi.atomOfOrb[m]) {
											size++;
										}
									}

								}
							}
						}
					}
					else if (mi.atomOfOrb[j] == mi.atomOfOrb[k]) {
						size++;

						for (int l : mi.missingOfAtom[mi.atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : mi.missingOfAtom[mi.atomOfOrb[j]]) {
									if (m > -1) {
										if (mi.atomOfOrb[l] == mi.atomOfOrb[m]) {
											size++;
										}
									}

								}
							}
						}
					}
					else {
						for (int l : mi.orbsOfAtom[mi.atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : mi.orbsOfAtom[mi.atomOfOrb[k]]) {
									if (m > -1) {
										size++;
									}
								}
							}
						}
					}
				}
			}

			return size;
		}

		@SuppressWarnings("DuplicatedCode")
		public static int findNCoulombInts(MoleculeInfo mi) {
			int size = 0;

			for (int j = 0; j < mi.nOrbitals; j++) {
				for (int k = j; k < mi.nOrbitals; k++) {
					if (j == k) {
						for (int l : mi.orbsOfAtom[mi.atomOfOrb[j]]) {
							if (l > -1) {
								size++;
							}
						}
						for (int l : mi.missingOfAtom[mi.atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : mi.missingOfAtom[mi.atomOfOrb[j]]) {
									if (m > -1) {
										if (mi.atomOfOrb[l] == mi.atomOfOrb[m]) {
											size++;
										}
									}
								}
							}
						}
					}
					else if (mi.atomOfOrb[j] == mi.atomOfOrb[k]) {
						size++;
						for (int l : mi.missingOfAtom[mi.atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : mi.missingOfAtom[mi.atomOfOrb[j]]) {
									if (m > -1) {
										if (mi.atomOfOrb[l] == mi.atomOfOrb[m]) {
											size++;
										}
									}

								}
							}
						}
					}
				}
			}

			return size;
		}

		public static int findNExchangeInts(MoleculeInfo mi) {
			int size = 0;

			for (int j = 0; j < mi.nOrbitals; j++) {
				for (int k = j; k < mi.nOrbitals; k++) {
					if (j == k) {
						for (int l : mi.orbsOfAtom[mi.atomOfOrb[j]]) {
							if (l > -1) {
								size++;
							}
						}
					}
					else if (mi.atomOfOrb[j] == mi.atomOfOrb[k]) {
						size++;
					}
					else {
						for (int l : mi.orbsOfAtom[mi.atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : mi.orbsOfAtom[mi.atomOfOrb[k]]) {
									if (m > -1) {
										size++;
									}
								}
							}
						}
					}
				}
			}

			return size;
		}

		/**
		 * npMap is the only model dependent field in RawMolecule, so it should be decided at build time.
		 *
		 * @param npMap "Map" of param numbers needed for differentiation. Key is Z.
		 * @return A complete, valid RawMolecule object.
		 */
		public RawMolecule build(int[][] npMap) throws IOException {
			if (atoms == null) throw new IllegalArgumentException("Atoms cannot be null!");
			if (datum == null) throw new IllegalArgumentException("Datum cannot be null!");

			RawMolecule rm = new RawMolecule();
			rm.atoms = atoms;
			rm.expGeom = expGeom;
			rm.datum = datum;
			rm.charge = charge;
			rm.mult = mult;
			rm.restricted = restricted;

			StringBuilder nameBuilder = new StringBuilder();
			Map<String, Integer> nameOccurrences = new TreeMap<>(Collections.reverseOrder());
			List<Integer> atomicNumbers = new ArrayList<>();

			// Map of unique Zs and their corresponding param numbers needed for differentiation.
			Map<Integer, int[]> uniqueNPs = new TreeMap<>();

			for (RawAtom a : rm.atoms) {
				atomicNumbers.add(a.Z);

				AtomProperties ap = AtomHandler.atoms[a.Z];
				rm.nElectrons += ap.getQ();
				rm.nOrbitals += ap.getOrbitals().length;

				if (!uniqueNPs.containsKey(ap.getZ())) {
					uniqueNPs.put(ap.getZ(), npMap[ap.getZ()]);
				}

				String name = ap.getName();
				if (!nameOccurrences.containsKey(name)) nameOccurrences.put(name, 1);
				else nameOccurrences.put(name, nameOccurrences.get(name) + 1);
			}

			for (String key : nameOccurrences.keySet()) nameBuilder.append(key).append(nameOccurrences.get(key));
			rm.name = nameBuilder + "_" + rm.charge + "_" + (rm.restricted ? "RHF" : "UHF");

			rm.nElectrons -= rm.charge;

			rm.atomicNumbers = Utils.toInts(atomicNumbers);

			rm.mats = Utils.toInts(uniqueNPs.keySet());

			rm.mnps = new int[uniqueNPs.size()][];
			int i = 0;
			for (int[] np : uniqueNPs.values()) {
				rm.mnps[i] = np;
				i++;
			}


			AtomHandler.populateAtoms();

			rm.orbsOfAtom = new int[rm.atoms.length][4];
			rm.atomOfOrb = new int[rm.nOrbitals];
			rm.missingOfAtom = new int[rm.atoms.length][4 * (rm.atoms.length - 1)];
			for (int[] index : rm.missingOfAtom) Arrays.fill(index, -1);

			int overallOrbitalIndex = 0;
			for (int atomIndex = 0, atomsLength = rm.atoms.length; atomIndex < atomsLength; atomIndex++) {
				AtomProperties atom = AtomHandler.atoms[rm.atoms[atomIndex].Z];
				int orbitalIndex = 0;
				for (OrbitalProperties ignored : atom.getOrbitals()) {
					rm.atomOfOrb[overallOrbitalIndex] = atomIndex;
					rm.orbsOfAtom[atomIndex][orbitalIndex] = overallOrbitalIndex;
					overallOrbitalIndex++;
					orbitalIndex++;
				}

				if (atom.getZ() == 1) {
					rm.orbsOfAtom[atomIndex][1] = -1;
					rm.orbsOfAtom[atomIndex][2] = -1;
					rm.orbsOfAtom[atomIndex][3] = -1;
				}
			}

			for (int j = 0; j < rm.atoms.length; j++) {
				int[] nums = new int[]{rm.orbsOfAtom[j][0], rm.orbsOfAtom[j][1],
						rm.orbsOfAtom[j][2], rm.orbsOfAtom[j][3]};
				int counter = 0;
				for (int k = 0; k < rm.nOrbitals; k++) {
					if (k != nums[0] && k != nums[1] && k != nums[2] && k != nums[3]) {
						rm.missingOfAtom[j][counter] = k;
						counter++;
					}
				}
			}

			if (rm.restricted) {
				rm.nIntegrals = findNIntegrals(rm);
			}
			else {
				rm.nCoulombInts = findNCoulombInts(rm);
				rm.nExchangeInts = findNExchangeInts(rm);
			}

			return rm;
		}

		public int getIndex() {
			return index;
		}

		public void setIndex(int index) {
			this.index = index;
		}

		public RawAtom[] getAtoms() {
			return atoms;
		}

		public void setAtoms(RawAtom[] atoms) {
			this.atoms = atoms;
		}

		public RawAtom[] getExpGeom() {
			return expGeom;
		}

		public void setExpGeom(RawAtom[] expGeom) {
			this.expGeom = expGeom;
		}

		public double[] getDatum() {
			return datum;
		}

		public void setDatum(double[] datum) {
			this.datum = datum;
		}

		public int getCharge() {
			return charge;
		}

		public void setCharge(int charge) {
			this.charge = charge;
		}

		public int getMult() {
			return mult;
		}

		public void setMult(int mult) {
			this.mult = mult;
		}

		public boolean isRestricted() {
			return restricted;
		}

		public void setRestricted(boolean restricted) {
			this.restricted = restricted;
		}

		public void setAtoms(int[] Zs, double[][] coords) {
			atoms = toRawAtoms(Zs, coords);
		}

		public void setExpGeom(int[] Zs, double[][] coords) {
			expGeom = toRawAtoms(Zs, coords);
		}

		private RawAtom[] toRawAtoms(int[] Zs, double[][] coords) {
			if (Zs.length != coords.length) {
				throw new IllegalArgumentException("Zs and coordinates length mismatch!");
			}

			int length = Zs.length;

			RawAtom[] ras = new RawAtom[length];

			for (int i = 0; i < length; i++) {
				ras[i] = new RawAtom();
				ras[i].Z = Zs[i];
				ras[i].coords = coords[i];
			}

			return ras;
		}
	}
}
