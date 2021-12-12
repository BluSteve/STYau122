package runcycle.input;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import nddo.am1.AM1Atom;
import nddo.mndo.MNDOAtom;
import scf.AtomHandler;
import scf.Model;
import tools.Utils;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

public class InputHandler {
	/**
	 * Converts geometry coordinates into bohr units before it's used by the
	 * program.
	 *
	 * @param filename filename without extension
	 * @return rawinput object in bohr units
	 */
	public static RawInput processInput(String filename) throws IOException {
		RawInput ri = new Gson().fromJson(new FileReader(filename + ".json"),
				RawInput.class);
		for (RawMolecule rm : ri.molecules) {
			for (int i = 0; i < rm.atoms.length; i++) {
				rm.atoms[i].coords = Utils.bohr(rm.atoms[i].coords);
				if (rm.expGeom != null)
					rm.expGeom[i].coords =
							Utils.bohr(rm.expGeom[i].coords);
			}
		}
		return ri;
	}

	/**
	 * Converts from bohr to armstrongs before outputting. Should always be
	 * used instead of outputJSON for any "processed" input.
	 *
	 * @param ri       rawinput object in bohr units
	 * @param filename filename without extension
	 */
	public static void outputInput(RawInput ri, String filename)
			throws IOException {
		for (RawMolecule rm : ri.molecules) {
			for (int i = 0; i < rm.atoms.length; i++) {
				rm.atoms[i].coords = Utils.debohr(rm.atoms[i].coords);
				if (rm.expGeom != null) rm.expGeom[i].coords =
						Utils.debohr(rm.expGeom[i].coords);
			}
		}
		outputJSON(ri, filename);
	}

	/**
	 * Outputs input json file after computing hash. File output is only done here.
	 *
	 * @param ri       the rawinput object to transcribe as is
	 * @param filename filename without extension
	 */
	private static void outputJSON(RawInput ri, String filename)
			throws IOException {
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.serializeNulls().create();
		String jsoned = gson.toJson(ri);
		String hash = Utils.getHash(jsoned);
		ri.hash = hash; // as such, the hash of the final file is different

		FileWriter fw = new FileWriter(filename + ".json");
		Files.createDirectories(Path.of("pastinputs"));
		FileWriter fwarchive = new FileWriter(
				"pastinputs/" + filename + "-" + hash + ".json");
		gson.toJson(ri, fwarchive);
		gson.toJson(ri, fw);
		fw.close();
		fwarchive.close();
	}

	/**
	 * Converts from .txt files to json
	 *
	 * @param filename filename with extension
	 */
	private static void convertFromTXT(String filename) throws IOException {
		AtomHandler.populateAtoms();
		RawInput ri = new RawInput();
		List<String> lines = Files.readAllLines(Path.of(filename));
		List<String> datums = Files.readAllLines(Path.of("reference.txt"));
		List<String> hud = Files.readAllLines(Path.of("mndohessian.txt"));
		String[] mndoParamsS = Files.readAllLines(Path.of("mndoparams.txt")).get(0).split(", ");
		double[] mp = Utils.toDoubles(mndoParamsS);

		ri.model = Model.MNDO; // TODO change this for AM1
		ri.trainingSet = lines.get(0).split("=")[1];
		ri.atomTypes = new int[ri.trainingSet.length()];
		for (int x = 0; x < ri.trainingSet.length(); x++) {
			ri.atomTypes[x] = AtomHandler.atomsMap
					.get(Character.toString(ri.trainingSet.charAt(x)))
					.getZ();
		}


		// Params
		int[][] neededParams = new int[ri.atomTypes.length][];
		int[][] npMap = new int[Utils.maxAtomNum][];
		int i = 0;
		for (int atomType : ri.atomTypes) {
			switch (ri.model) {
				case MNDO:
					if (atomType == 1) neededParams[i] = npMap[atomType] = MNDOAtom.T1ParamNums;
					else neededParams[i] = npMap[atomType] = MNDOAtom.T2ParamNums;
					break;
				case AM1:
					if (atomType == 1) neededParams[i] = npMap[atomType] = AM1Atom.HParamNums;
					if (atomType == 5) neededParams[i] = npMap[atomType] = AM1Atom.NParamNums;
					if (atomType == 6) neededParams[i] = npMap[atomType] = AM1Atom.CParamNums;
					if (atomType == 8) neededParams[i] = npMap[atomType] = AM1Atom.OParamNums;
					break;
			}
			i++;
		}
		ri.neededParams = neededParams;

		int PARAMLENGTH;
		switch (ri.model) {
			case MNDO:
				PARAMLENGTH = 13;
				break;
			case AM1:
				PARAMLENGTH = 25;
				break;
			default:
				throw new IllegalStateException("Unexpected value: " + ri.model);
		}

		ri.params = new RawParams();
		ri.params.nddoParams = new double[ri.atomTypes.length][];
		ri.params.lastHessian = Utils.toDoubles(hud.get(0).split(":")[1].strip().split(","));
		ri.params.lastGradient = Utils.toDoubles(hud.get(1).split(":")[1].strip().split(","));
		ri.params.lastDir = Utils.toDoubles(hud.get(2).split(":")[1].strip().split(","));
		for (int x = 0; x < ri.atomTypes.length; x++)
			ri.params.nddoParams[x] = Arrays.copyOfRange(mp, x * PARAMLENGTH, x * PARAMLENGTH + PARAMLENGTH);


		// Molecules
		ArrayList<RawMolecule> moleculesL = new ArrayList<>();
		i = 1;
		int datumi = 0;
		int moleculei = 0;

		while (i < lines.size()) {
			RawMolecule.RMBuilder builder = new RawMolecule.RMBuilder();

			builder.setRestricted(lines.get(i).equals("RHF"));
			i++;

			builder.setCharge(Integer.parseInt(lines.get(i).split("=")[1]));
			i++;

			builder.setMult(Integer.parseInt(lines.get(i).split("=")[1]));
			i++;

			builder.setIndex(moleculei);

			// Atoms
			ArrayList<RawAtom> atomsL = new ArrayList<>();

			while (!lines.get(i).equals("---") && !lines.get(i).equals("EXPGEOM")) {
				addAtom(atomsL, lines, i);
				i++;
			}
			RawAtom[] atoms = new RawAtom[atomsL.size()];
			for (int p = 0; p < atomsL.size(); p++) atoms[p] = atomsL.get(p);
			builder.setAtoms(atoms);


			// expGeom
			ArrayList<RawAtom> expGeomL = new ArrayList<>();
			RawAtom[] expGeom = null;
			if (lines.get(i).equals("EXPGEOM")) {
				i++;
				while (!lines.get(i).equals("---")) {
					addAtom(expGeomL, lines, i);
					i++;
				}
				expGeom = new RawAtom[expGeomL.size()];
				for (int p = 0; p < expGeomL.size(); p++)
					expGeom[p] = expGeomL.get(p);
			}
			builder.setExpGeom(expGeom);

			if (expGeomL.size() != 0 && atomsL.size() != expGeomL.size())
				throw new IllegalArgumentException("Atom and expGeom size mismatch!");


			// Datum
			double[] datum = new double[3];

			datum[0] = Double.parseDouble(datums.get(datumi).split(" ")[1]);
			datumi++;

			String[] ss = datums.get(datumi).split(" ");
			if (ss.length > 1) datum[1] = Double.parseDouble(ss[1]);
			datumi++;

			ss = datums.get(datumi).split(" ");
			if (ss.length > 1) datum[2] = Double.parseDouble(ss[1]);
			datumi += 2;

			builder.setDatum(datum);


			moleculesL.add(builder.build(npMap));
			i++;
		}

		ri.molecules = new RawMolecule[moleculesL.size()];
		for (int p = 0; p < moleculesL.size(); p++) {
			moleculesL.get(p).index = p;
			ri.molecules[p] = moleculesL.get(p);
		}
		ri.nMolecules = ri.molecules.length;

		outputJSON(ri, "input");
	}

	private static void outputSubset(int... ids) throws IOException {
		RawInput ri = processInput("input");
		ri.nMolecules = ids.length;
		RawMolecule[] newrms = new RawMolecule[ids.length];
		int ni = 0;
		for (int i = 0; i < ri.molecules.length; i++) {
			for (int id : ids) {
				if (ri.molecules[i].index == id) {
					newrms[ni] = ri.molecules[i];
					ni++;
					break;
				}
			}
		}
		ri.molecules = newrms;
		outputInput(ri, "subset");
	}

	private static void addAtom(List<RawAtom> rawAtoms, List<String> lines, int i) {
		RawAtom ra = new RawAtom();
		StringTokenizer t = new StringTokenizer(lines.get(i), " ");
		t.nextToken();
		String name = t.nextToken();
		ra.Z = AtomHandler.atomsMap.get(name).getZ();
		ra.coords = new double[3];
		for (int q = 0; q < 3; q++)
			ra.coords[q] = Double.parseDouble(t.nextToken());
		rawAtoms.add(ra);
	}

	public static void main(String[] args) throws IOException {
		convertFromTXT("inputtesting.txt");
		outputSubset(9);
	}
}
