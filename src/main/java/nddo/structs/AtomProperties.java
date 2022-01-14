package nddo.structs;

import com.google.gson.Gson;
import nddo.Constants;
import tools.Utils;

import java.util.Arrays;
import java.util.HashMap;

public class AtomProperties { // Only 119 of these, immutable
	private int Z, Q;
	private String name;
	private double mass, heat;	private static AtomProperties[] atoms = getAtoms();
	private OrbitalProperties[] orbitals;

	private static void populateAtoms() {
		atoms = new AtomProperties[Constants.MAX_ATOM_NUM];
		atomsMap = new HashMap<>();

		try {
			String json = Utils.getResource("atom-properties.json");

			Gson gson = new Gson();

			AtomProperties[] unindexedAtoms = gson.fromJson(json, AtomProperties[].class);

			for (AtomProperties atom : unindexedAtoms) {
				atom.setOrbitals(OrbitalProperties.generateOrbitals(atom.getPeriod()));

				getAtoms()[atom.getZ()] = atom;
				getAtomsMap().put(atom.getName(), atom);
			}
		} catch (NullPointerException e) {
			throw new RuntimeException("atom-properties.json not found!");
		}
	}

	public static AtomProperties[] getAtoms() {
		if (atoms == null) populateAtoms();

		return atoms;
	}

	public static HashMap<String, AtomProperties> getAtomsMap() {
		if (atomsMap == null) populateAtoms();

		return atomsMap;
	}	private static HashMap<String, AtomProperties> atomsMap = getAtomsMap();

	public int getPeriod() {
		int period = 1;
		if (Z > 85) period = 7;
		else if (Z > 53) period = 6;
		else if (Z > 35) period = 5;
		else if (Z > 17) period = 4;
		else if (Z > 9) period = 3;
		else if (Z > 1) period = 2;
		return period;
	}

	public int getZ() {
		return Z;
	}

	public int getQ() {
		return Q;
	}

	public double getHeat() {
		return heat * Constants.HEATCONV;
	}

	public double getMass() {
		return mass;
	}

	public String getName() {
		return name;
	}

	public OrbitalProperties[] getOrbitals() {
		return orbitals;
	}

	void setOrbitals(OrbitalProperties[] orbitals) {
		this.orbitals = orbitals;
	}

	@Override
	public String toString() {
		return "AtomProperties{" +
				", Z=" + Z +
				", Q=" + Q +
				", name='" + name + '\'' +
				", mass=" + mass +
				", heat=" + heat +
				", orbitals=" + Arrays.toString(orbitals) +
				'}';
	}






}
