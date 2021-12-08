package runcycle.input;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import nddo.solution.Solution;
import nddo.am1.AM1Params;
import nddo.mndo.MNDOParams;
import scf.AtomHandler;
import scf.AtomProperties;
import tools.Utils;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

import static nddo.mndo.MNDOParams.T1ParamNums;
import static nddo.mndo.MNDOParams.T2ParamNums;

public class InputHandler {
	/**
	 * Converts geometry coordinates into bohr units before it's used by the
	 * program.
	 *
	 * @param filename filename without extension
	 * @return rawinput object in bohr units
	 */
	public static RawInput processInput(String filename) throws IOException {
		RawInput ri = (new Gson()).fromJson(new FileReader(filename + ".json"),
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
		List<String> hud = Files.readAllLines(Path.of(
				"mndohessian.txt"));
		String[] mndoParamsS =
				Files.readAllLines(Path.of("mndoparams.txt")).get(0)
						.split(", ");
		double[] mp = Utils.toDoubles(mndoParamsS);

		ri.model = "mndo"; // TODO change this for AM1
		ri.trainingSet = lines.get(0).split("=")[1];
		ri.atomTypes = new int[ri.trainingSet.length()];
		for (int x = 0; x < ri.trainingSet.length(); x++) {
			ri.atomTypes[x] = AtomHandler.atomsMap
					.get(Character.toString(ri.trainingSet.charAt(x)))
					.getZ();
		}
		int[][] neededParams = new int[ri.atomTypes.length][];
		int w = 0;
		for (int atomType : ri.atomTypes) {
			switch (ri.model) {
				case "mndo":
					if (atomType == 1)
						neededParams[w] = MNDOParams.T1ParamNums;
					else neededParams[w] = MNDOParams.T2ParamNums;
					break;
				case "am1":
					if (atomType == 1)
						neededParams[w] = AM1Params.HParamNums;
					if (atomType == 5) neededParams[w] =
							AM1Params.NParamNums;
					if (atomType == 6) neededParams[w] =
							AM1Params.CParamNums;
					if (atomType == 8) neededParams[w] =
							AM1Params.OParamNums;
					break;
			}
			w++;
		}
		ri.neededParams = neededParams;
		int PARAMLENGTH = 13;

		ri.params = new RawParams();
		ri.params.nddoParams = new double[ri.atomTypes.length][];
		ri.params.lastHessian =
				Utils.toDoubles(
						hud.get(0).split(":")[1].strip().split(","));
		ri.params.lastGradient =
				Utils.toDoubles(
						hud.get(1).split(":")[1].strip().split(","));
		ri.params.lastDir =
				Utils.toDoubles(
						hud.get(2).split(":")[1].strip().split(","));
		for (int x = 0; x < ri.atomTypes.length; x++) {
			ri.params.nddoParams[x] =
					Arrays.copyOfRange(mp, x * PARAMLENGTH,
							x * PARAMLENGTH +
									PARAMLENGTH);
		}

		ArrayList<RawMolecule> moleculesL = new ArrayList<>();

		int i = 1;
		try {
			while (i < lines.size()) {
				RawMolecule rm = new RawMolecule();
				ArrayList<RawAtom> atomsL = new ArrayList<>();
				ArrayList<RawAtom> expGeomL = new ArrayList<>();

				rm.restricted = lines.get(i).equals("RHF");
				i++;

				rm.charge = Integer.parseInt(lines.get(i).split("=")[1]);
				i++;

				rm.mult = Integer.parseInt(lines.get(i).split("=")[1]);
				i++;

				while (!lines.get(i).equals("---") &&
						!lines.get(i).equals("EXPGEOM")) {
					addAtom(atomsL, lines, i);
					i++;
				}

				StringBuilder nameBuilder = new StringBuilder();
				TreeMap<String, Integer> nameOccurrences =
						new TreeMap<>(Collections.reverseOrder());
				ArrayList<Integer> tempZs =
						new ArrayList<>(Utils.maxAtomNum);
				ArrayList<Integer> atomicNumbers = new ArrayList<>();
				HashMap<Integer, int[]> tempNPs = new HashMap<>();
				for (RawAtom a : atomsL) {
					atomicNumbers.add(a.Z);
					AtomProperties ap = AtomHandler.atomsMap.get(a.name);
					rm.nElectrons += ap.getQ();
					rm.nOrbitals += ap.getOrbitals().length;
					if (!tempZs.contains(ap.getZ())) {
						tempZs.add(ap.getZ());
						if (ri.model.equals("mndo")) {
							if (ap.getZ() == 1)
								tempNPs.put(ap.getZ(), T1ParamNums);
							else tempNPs.put(ap.getZ(), T2ParamNums);
						}
						// todo implement for am1
					}
					if (!nameOccurrences.containsKey(a.name))
						nameOccurrences.put(a.name, 1);
					else
						nameOccurrences.put(a.name,
								nameOccurrences.get(a.name) + 1);
				}
				rm.nElectrons -= rm.charge;
				rm.atomicNumbers = Utils.toInts(atomicNumbers);
				for (String key : nameOccurrences.keySet()) {
					nameBuilder.append(key)
							.append(nameOccurrences.get(key));
				}
				Collections.sort(tempZs);
				int[] moleculeATs = Utils.toInts(tempZs);
				int[][] moleculeNPs = new int[tempNPs.size()][];
				for (int j = 0; j < moleculeNPs.length; j++) {
					moleculeNPs[j] = tempNPs.get(tempZs.get(j));
				}
				rm.mnps = moleculeNPs;
				rm.mats = moleculeATs;
				String ruhf = rm.restricted ? "RHF" : "UHF";
				rm.name = nameBuilder + "_" + rm.charge + "_" + ruhf;

				if (lines.get(i).equals("EXPGEOM")) {
					i++;
					while (!lines.get(i).equals("---")) {
						addAtom(expGeomL, lines, i);
						i++;
					}
					rm.expGeom = new RawAtom[atomsL.size()];
					for (int p = 0; p < expGeomL.size(); p++)
						rm.expGeom[p] = expGeomL.get(p);
				}
				else {
					rm.expGeom = null;
				}

				rm.atoms = new RawAtom[atomsL.size()];
				for (int p = 0; p < atomsL.size(); p++)
					rm.atoms[p] = atomsL.get(p);

				int[] nIntegrals = Solution.getNIntegrals(rm);
				if (rm.restricted)
					rm.nIntegrals = nIntegrals[0];
				else {
					rm.nCoulombInts = nIntegrals[0];
					rm.nExchangeInts = nIntegrals[1];
				}
				moleculesL.add(rm);
				i++;
			}
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		}
		try {
			i = 0;
			int j = 0;
			while (i < datums.size() && j < moleculesL.size()) {
				String[] ss;
				double[] datum = new double[3];
				datum[0] = Double.parseDouble(datums.get(i).split(" ")[1]);
				i++;
				ss = datums.get(i).split(" ");
				if (ss.length == 2) datum[1] = Double.parseDouble(ss[1]);
				i++;
				ss = datums.get(i).split(" ");
				if (ss.length == 2) datum[2] = Double.parseDouble(ss[1]);
				i += 2;
				moleculesL.get(j).datum = datum;
				j++;
			}
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		}
		ri.molecules = new RawMolecule[moleculesL.size()];
		for (int p = 0; p < moleculesL.size(); p++)
			ri.molecules[p] = moleculesL.get(p);
		for (int j = 0; j < ri.molecules.length; j++)
			ri.molecules[j].index = j;
		ri.nMolecules = ri.molecules.length;

		outputJSON(ri, "input");
	}

	private static void outputSubset(int... ids)
			throws IOException {
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

	private static void addAtom(List<RawAtom> rawAtoms, List<String> lines,
								int i) {
		RawAtom ra = new RawAtom();
		StringTokenizer t = new StringTokenizer(lines.get(i), " ");
		t.nextToken();
		ra.name = t.nextToken();
		ra.Z = AtomHandler.atomsMap.get(ra.name).getZ();
		for (int q = 0; q < 3; q++)
			ra.coords[q] = Double.parseDouble(t.nextToken());
		rawAtoms.add(ra);
	}

	public static void main(String[] args) throws IOException {
		convertFromTXT("inputtesting.txt");
//		outputSubset(new int[]{234, 237, 258, 260, 261, 270, 276, 285});
		outputSubset(9);
	}
}
