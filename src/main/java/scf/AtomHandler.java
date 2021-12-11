package scf;


import com.google.gson.Gson;
import tools.Utils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;

public class AtomHandler {
	public static AtomProperties[] atoms = new AtomProperties[Utils.maxAtomNum];
	public static HashMap<String, AtomProperties> atomsMap = new HashMap<>();

	public static void populateAtoms() throws IOException {
		String json = Files.readString(Path.of("atom-properties.json"));
		Gson gson = new Gson();

		AtomProperties[] unindexedAtoms = gson.fromJson(json, AtomProperties[].class);

		for (AtomProperties atom : unindexedAtoms) {
			atom.setOrbitals(OrbitalProperties.generateOrbitals(atom.getPeriod()));

			atoms[atom.getZ()] = atom;
			atomsMap.put(atom.getName(), atom);
		}
	}
}
