package runcycle.input;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import nddoparam.Solution;
import nddoparam.am1.AM1Params;
import nddoparam.mndo.MNDOParams;
import scf.AtomHandler;
import scf.AtomProperties;
import tools.Utils;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

import static nddoparam.mndo.MNDOParams.T1ParamNums;
import static nddoparam.mndo.MNDOParams.T2ParamNums;

public class InputHandler {
	public static RawInput ri;

	public static void processInput(String path) {
		try {
			ri = (new Gson()).fromJson(new FileReader(path + ".json"),
					RawInput.class);
			for (RawMolecule rm : ri.molecules) {
				for (int i = 0; i < rm.atoms.length; i++) {
					rm.atoms[i].coords = Utils.bohr(rm.atoms[i].coords);
					if (rm.expGeom != null)
						rm.expGeom[i].coords =
								Utils.bohr(rm.expGeom[i].coords);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void updateInput(RawInput ri, String input) {
		for (RawMolecule rm : ri.molecules) {
			for (int i = 0; i < rm.atoms.length; i++) {
				rm.atoms[i].coords = Utils.debohr(rm.atoms[i].coords);
				if (rm.expGeom != null) rm.expGeom[i].coords =
						Utils.debohr(rm.expGeom[i].coords);
			}
		}
		makeInput(ri, input);
	}

	public static void makeInput(RawInput ri, String input) {
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.serializeNulls().create();
		String jsoned = gson.toJson(ri);
		String hash = Utils.getHash(jsoned);
		ri.hash = hash; // as such, the hash of the final file is different

		try {
			FileWriter fw = new FileWriter(input + ".json");
			Files.createDirectories(Path.of("pastinputs"));
			FileWriter fwarchive = new FileWriter(
					"pastinputs/" + input + "-" + hash + ".json");
			gson.toJson(ri, fwarchive);
			gson.toJson(ri, fw);
			fw.close();
			fwarchive.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static void convertFromTXT(String input) {
		AtomHandler.populateAtoms();
		ri = new RawInput();
		try {
			List<String> lines = Files.readAllLines(Path.of(input));
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
				while (i < datums.size()) {
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

			makeInput(ri, "input");
		} catch (
				IOException e) {
			e.printStackTrace();
		}

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

	public static void main(String[] args) {
		InputHandler.convertFromTXT("inputtesting.txt");
	}
}
