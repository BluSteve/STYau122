package nddo.structs;

import com.google.gson.Gson;
import nddo.Constants;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;

public class AtomProperties { // Only 119 of these, immutable

	private static AtomProperties[] atoms = getAtoms();
	private static HashMap<String, AtomProperties> atomsMap = getAtomsMap();

	private static void populateAtoms() {
		atoms = new AtomProperties[Constants.maxAtomNum];
		atomsMap = new HashMap<>();

		try {
			Gson gson = new Gson();
			String json = Files.readString(Path.of("atom-properties.json"));

			AtomProperties[] unindexedAtoms = gson.fromJson(json, AtomProperties[].class);

			for (AtomProperties atom : unindexedAtoms) {
				atom.setOrbitals(OrbitalProperties.generateOrbitals(atom.getPeriod()));

				getAtoms()[atom.getZ()] = atom;
				getAtomsMap().put(atom.getName(), atom);
			}
		} catch (IOException e) {
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
	}
	private static final double HEATCONV = 4.3363E-2;
	private int Z, Q;
	private String name;
	private double mass, heat;
	private OrbitalProperties[] orbitals;

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
		return heat * HEATCONV;
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
