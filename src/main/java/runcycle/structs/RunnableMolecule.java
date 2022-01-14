package runcycle.structs;

import nddo.NDDOAtom;
import nddo.NDDOParams;
import nddo.solution.Solution;
import nddo.structs.AtomProperties;
import nddo.structs.MoleculeInfo;
import org.ejml.simple.SimpleMatrix;
import runcycle.State;
import tools.Utils;

import java.util.*;

public final class RunnableMolecule extends MoleculeInfo { // mid-level runnable molecule representation
	public final Atom[] atoms, expGeom;
	public final double[] datum;

	// MoleculeInfo will ever change, even across runs!
	public RunnableMolecule(MoleculeInfo mi, Atom[] atoms, Atom[] expGeom, double[] datum) {
		super(mi);

		if (atoms == null) throw new NullPointerException("Atoms cannot be null!");
		if (datum == null) throw new NullPointerException("Datum cannot be null!");
		if (expGeom != null && atoms.length != expGeom.length)
			throw new IllegalArgumentException("Atom and expGeom size mismatch!");

		this.atoms = atoms;
		this.expGeom = expGeom;
		this.datum = datum;
	}

	@Override
	public String toString() {
		return Serializer.gson.toJson(this);
	}

	public static class RMBuilder {
		public int index;
		public String name;
		public boolean restricted = true;
		public Boolean useEdiis;
		public double[][] densityMatrices, densityMatricesExp;
		public int charge, mult;
		public double[] datum;
		public Atom[] atoms, expGeom;

		public RunnableMolecule build(int[] atomTypes, int[][] neededParams,
									  NDDOParams[] npMap) { // can use neededParams map instead
			Arrays.sort(atoms, Comparator.comparingInt(o -> -o.Z));
			if (expGeom != null) Arrays.sort(expGeom, Comparator.comparingInt(o -> -o.Z));

			MoleculeInfo.MIBuilder miBuilder = new MoleculeInfo.MIBuilder();

			miBuilder.index = index;
			miBuilder.charge = charge;
			miBuilder.mult = mult;
			miBuilder.restricted = restricted;

			StringBuilder nameBuilder = new StringBuilder();
			Map<String, Integer> nameOccurrences = new TreeMap<>(
					Comparator.comparingInt(o -> AtomProperties.getAtomsMap().get(o).getZ()));

			List<Integer> atomicNumbers = new ArrayList<>();

			// Map of unique Zs and their corresponding param numbers needed for differentiation.
			Map<Integer, int[]> uniqueNPs = new TreeMap<>();

			for (Atom a : atoms) {
				atomicNumbers.add(a.Z);

				AtomProperties ap = AtomProperties.getAtoms()[a.Z];
				miBuilder.nElectrons += ap.getQ();
				miBuilder.nOrbitals += ap.getOrbitals().length;

				if (!uniqueNPs.containsKey(ap.getZ())) {
					// searches through until it finds corresponding neededParams
					for (int i = 0; i < atomTypes.length; i++) {
						if (atomTypes[i] == ap.getZ()) {
							uniqueNPs.put(ap.getZ(), neededParams[i]);
							break;
						}
					}
				}

				if (name == null) {
					String atomName = ap.getName();
					if (!nameOccurrences.containsKey(atomName)) nameOccurrences.put(atomName, 1);
					else nameOccurrences.put(atomName, nameOccurrences.get(atomName) + 1);
				}
			}

			if (name == null) { // default name
				for (String key : nameOccurrences.keySet()) nameBuilder.append(key).append(nameOccurrences.get(key));
				miBuilder.name = nameBuilder.toString();
			}
			else miBuilder.name = name;

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
				AtomProperties atom = AtomProperties.getAtoms()[atoms[atomIndex].Z];

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

			MoleculeInfo temp = miBuilder.build();

			miBuilder.densityMatrices = densityMatrices;
			miBuilder.densityMatricesExp = densityMatricesExp;

			if (npMap != null && (useEdiis == null || densityMatrices == null || densityMatricesExp == null)) {
				NDDOAtom[] nddoAtoms = State.getConverter().convert(atoms, npMap);

				miBuilder.useEdiis = useEdiis == null ? Solution.shouldEdiis(temp, nddoAtoms) : useEdiis;

				if (densityMatrices == null && densityMatricesExp == null) {
					temp = miBuilder.build();
					miBuilder.densityMatrices = getDoubleDM(temp, nddoAtoms);
					if (expGeom != null) {
						NDDOAtom[] expAtoms = State.getConverter().convert(expGeom, npMap);
						miBuilder.densityMatricesExp = getDoubleDM(temp, expAtoms);
					}
				}
			}

			return new RunnableMolecule(miBuilder.build(), atoms, expGeom, datum);
		}

		private double[][] getDoubleDM(MoleculeInfo temp, NDDOAtom[] nddoAtoms) {
			SimpleMatrix[] dms = Solution.findDensityMatrices(temp, nddoAtoms);
			double[][] doubleDM = new double[dms.length][];
			for (int j = 0; j < dms.length; j++) {
				doubleDM[j] = dms[j].getDDRM().data;
			}

			return doubleDM;
		}

		public void setAtoms(int[] Zs, double[][] coords) {
			atoms = toAtoms(Zs, coords);
		}

		public void setExpGeom(int[] Zs, double[][] coords) {
			expGeom = toAtoms(Zs, coords);
		}

		private Atom[] toAtoms(int[] Zs, double[][] coords) {
			if (Zs.length != coords.length) {
				throw new IllegalArgumentException("Zs and coordinates length mismatch!");
			}

			int length = Zs.length;

			Atom[] ras = new Atom[length];

			for (int i = 0; i < length; i++) {
				ras[i] = new Atom(Zs[i], coords[i]);
			}

			return ras;
		}
	}
}
