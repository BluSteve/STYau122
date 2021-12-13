package runcycle.input;

import nddo.structs.AtomHandler;
import nddo.structs.AtomProperties;
import nddo.structs.MoleculeInfo;
import tools.Utils;

import java.util.*;

public class RawMolecule extends MoleculeInfo { // mid-level runnable molecule representation
	public final RawAtom[] atoms, expGeom;
	public final double[] datum;

	private RawMolecule(MoleculeInfo mi, RawAtom[] atoms, RawAtom[] expGeom, double[] datum) {
		super(mi);
		this.atoms = atoms;
		this.expGeom = expGeom;
		this.datum = datum;
	}

	public static class RMBuilder {
		public int index;
		public boolean restricted = true;
		public int charge, mult;
		public double[] datum;
		public RawAtom[] atoms, expGeom;

		/**
		 * npMap is the only model dependent field in RawMolecule, so it should be decided at build time.
		 *
		 * @param npMap "Map" of param numbers needed for differentiation. Key is Z.
		 * @return A complete, valid RawMolecule object.
		 */
		public RawMolecule build(int[][] npMap) {
			if (atoms == null) throw new IllegalArgumentException("Atoms cannot be null!");
			if (datum == null) throw new IllegalArgumentException("Datum cannot be null!");

			MoleculeInfo.MIBuilder miBuilder = new MoleculeInfo.MIBuilder();

			miBuilder.index = index;
			miBuilder.charge = charge;
			miBuilder.mult = mult;
			miBuilder.restricted = restricted;

			StringBuilder nameBuilder = new StringBuilder();
			Map<String, Integer> nameOccurrences = new TreeMap<>(Collections.reverseOrder());
			List<Integer> atomicNumbers = new ArrayList<>();

			// Map of unique Zs and their corresponding param numbers needed for differentiation.
			Map<Integer, int[]> uniqueNPs = new TreeMap<>();

			for (RawAtom a : atoms) {
				atomicNumbers.add(a.Z);

				AtomProperties ap = AtomHandler.getAtoms()[a.Z];
				miBuilder.nElectrons += ap.getQ();
				miBuilder.nOrbitals += ap.getOrbitals().length;

				if (!uniqueNPs.containsKey(ap.getZ())) {
					uniqueNPs.put(ap.getZ(), npMap[ap.getZ()]);
				}

				String name = ap.getName();
				if (!nameOccurrences.containsKey(name)) nameOccurrences.put(name, 1);
				else nameOccurrences.put(name, nameOccurrences.get(name) + 1);
			}

			for (String key : nameOccurrences.keySet()) nameBuilder.append(key).append(nameOccurrences.get(key));
			miBuilder.name = nameBuilder + "_" + miBuilder.charge + "_" + (miBuilder.restricted ? "RHF" : "UHF");

			miBuilder.nElectrons -= miBuilder.charge;

			miBuilder.atomicNumbers = Utils.toInts(atomicNumbers);

			miBuilder.mats = Utils.toInts(uniqueNPs.keySet());

			miBuilder.mnps = new int[uniqueNPs.size()][];
			int i = 0;
			for (int[] np : uniqueNPs.values()) {
				miBuilder.mnps[i] = np;
				i++;
			}


			// Computes cached Solution info
			miBuilder.orbsOfAtom = new int[atoms.length][];
			miBuilder.atomOfOrb = new int[miBuilder.nOrbitals];
			miBuilder.missingOfAtom = new int[atoms.length][];

			int overallOrbitalIndex = 0;

			for (int atomIndex = 0; atomIndex < atoms.length; atomIndex++) {
				AtomProperties atom = AtomHandler.getAtoms()[atoms[atomIndex].Z];

				int olength = atom.getOrbitals().length;
				miBuilder.orbsOfAtom[atomIndex] = new int[olength];

				for (int orbitalIndex = 0; orbitalIndex < olength; orbitalIndex++) {
					miBuilder.atomOfOrb[overallOrbitalIndex] = atomIndex;

					miBuilder.orbsOfAtom[atomIndex][orbitalIndex] = overallOrbitalIndex;

					overallOrbitalIndex++;
				}
			}

			for (int j = 0; j < atoms.length; j++) {
				int[] missing = new int[miBuilder.nOrbitals - miBuilder.orbsOfAtom[j].length];

				int counter = 0;

				for (int k = 0; k < miBuilder.nOrbitals; k++) {
					// if k not in orbsOfAtom
					boolean in = false;
					for (int orb : miBuilder.orbsOfAtom[j]) {
						if (k == orb) {
							in = true;
							break;
						}
					}

					if (!in) {
						missing[counter] = k;
						counter++;
					}
				}

				miBuilder.missingOfAtom[j] = missing;
			}

			return new RawMolecule(miBuilder.build(), atoms, expGeom, datum);
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
