package scf;


import com.google.gson.Gson;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;

public class AtomHandler {
	public static AtomProperties[] atoms = new AtomProperties[119];
	public static HashMap<String, AtomProperties> atomsMap = new HashMap<>();

	public static void populateAtoms() throws IOException {
		String json = Files.readString(Path.of("atom-properties.json"));
		Gson gson = new Gson();

		AtomProperties[] unindexedAtoms = gson.fromJson(json, AtomProperties[].class);

		for (int i = 0; i < unindexedAtoms.length; i++) {
			AtomProperties atom = unindexedAtoms[i];

			atom.setIndex(i);
			atom.setOrbitals(OrbitalProperties.generateOrbitals(atom.getPeriod()));

			atoms[atom.getZ()] = atom;
			atomsMap.put(atom.getName(), atom);
		}
	}
}
