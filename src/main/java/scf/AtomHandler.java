package scf;


import com.google.gson.Gson;
import tools.Utils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;

public class AtomHandler {
	private static AtomProperties[] atoms;
	private static HashMap<String, AtomProperties> atomsMap;

	private static void populateAtoms() {
		atoms = new AtomProperties[Utils.maxAtomNum];
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
}
