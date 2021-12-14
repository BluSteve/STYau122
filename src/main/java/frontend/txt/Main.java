package frontend.txt;

import examples.am1.AM1Atom;
import frontend.json.RawInput;
import nddo.Constants;
import nddo.NDDOParams;
import nddo.mndo.MNDOAtom;
import nddo.structs.AtomProperties;
import runcycle.structs.Atom;
import runcycle.structs.Params;
import runcycle.structs.RunnableMolecule;
import tools.Utils;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

public class Main {
	public static void txtToText() throws IOException {
		List<String> lines = Files.readAllLines(Path.of("input.txt"));
		List<String> datums = Files.readAllLines(Path.of("reference.txt"));
		PrintWriter pw = new PrintWriter("molecules.txt");

		// Molecules
		ArrayList<RunnableMolecule> moleculesL = new ArrayList<>();
		int inputi = 1;
		int datumi = 0;

		while (inputi < lines.size()) {
			String rhf = lines.get(inputi);
			inputi++;

			String charge = lines.get(inputi);
			inputi++;

			String mult = lines.get(inputi);
			inputi++;

			// Atoms
			ArrayList<Atom> atomsL = new ArrayList<>();

			String atoms = "";
			while (!lines.get(inputi).equals("---") && !lines.get(inputi).equals("EXPGEOM")) {
				addAtom(atomsL, lines.get(inputi));
				atoms += lines.get(inputi) + "\n";
				inputi++;
			}

			Map<String, Integer> nameOccurrences = new TreeMap<>(
					Comparator.comparingInt(o -> AtomProperties.getAtomsMap().get(o).getZ()));

			for (Atom atom : atomsL) {
				AtomProperties ap = AtomProperties.getAtoms()[atom.Z];
				String name = ap.getName();

				if (!nameOccurrences.containsKey(name)) nameOccurrences.put(name, 1);
				else nameOccurrences.put(name, nameOccurrences.get(name) + 1);
			}

			String name = "";
			for (String key : nameOccurrences.keySet()) name += key + nameOccurrences.get(key);

			// expGeom
			String expGeom = "";
			if (lines.get(inputi).equals("EXPGEOM")) {
				while (!lines.get(inputi).equals("---")) {
					expGeom += lines.get(inputi) + "\n";
					inputi++;
				}
			}

			// Datum
			String[] datum = new String[3];

			datum[0] = "HF=" + datums.get(datumi).split(" ")[1];
			datumi++;

			String[] ss = datums.get(datumi).split(" ");
			if (ss.length > 1) datum[1] = "DIPOLE=" + ss[1];
			datumi++;

			ss = datums.get(datumi).split(" ");
			if (ss.length > 1) datum[2] = "IE=" + ss[1];
			datumi += 2;

			List<String> datumL = new ArrayList<>();
			for (String s : datum) {
				if (s != null) datumL.add(s);
			}


			pw.write(String.join(", ", name, rhf, charge, mult) + "\n");
			pw.write(String.join(", ", datumL) + "\n");
			pw.write(atoms);
			pw.write(expGeom);
			pw.write("---\n");


			inputi++;
		}

		pw.close();
	}

	private static void addAtom(List<Atom> Atoms, String line) {
		StringTokenizer t = new StringTokenizer(line, " ");
		t.nextToken();
		String name = t.nextToken();
		double[] coords = new double[3];
		for (int q = 0; q < 3; q++) coords[q] = Double.parseDouble(t.nextToken());
		Atom a = new Atom(AtomProperties.getAtomsMap().get(name).getZ(), coords);
		Atoms.add(a);
	}

	private static String[] splitCsvLine(String s) {
		return s.split(", *");
	}


	public static void convertFromTXT() throws IOException {
		List<String> pcsv = Files.readAllLines(Path.of("params.csv"));
		NDDOParams[] npMap = new NDDOParams[Constants.maxAtomNum];
		for (String line : pcsv) {
			String[] linea = splitCsvLine(line);
			npMap[AtomProperties.getAtomsMap().get(linea[0]).getZ()] =
					new NDDOParams(Utils.toDoubles(Arrays.copyOfRange(linea, 1, linea.length)));
		}

		List<String> pncsv = Files.readAllLines(Path.of("param-numbers.csv"));
		int[] atomTypes = new int[pncsv.size()];
		int[][] neededParams = new int[pncsv.size()][];
		int totalParamLength = 0;
		for (int i = 0; i < pncsv.size(); i++) {
			String[] linea = splitCsvLine(pncsv.get(i));
			atomTypes[i] = AtomProperties.getAtomsMap().get(linea[0]).getZ();
			neededParams[i] = Utils.toInts(Arrays.copyOfRange(linea, 1, linea.length));
			totalParamLength += linea.length - 1;
		}

		Params p = new Params(npMap, totalParamLength);

		Scanner s = new Scanner(new FileInputStream("molecules.txt"));

		// Molecules
		ArrayList<RunnableMolecule> moleculesL = new ArrayList<>();
		i = 1;
		int datumi = 0;
		int moleculei = 0;

		while (s.hasNextLine()) {
			RunnableMolecule.RMBuilder builder = new RunnableMolecule.RMBuilder();

			builder.restricted = params.get(i).equals("RHF");
			i++;

			builder.charge = Integer.parseInt(params.get(i).split("=")[1]);
			i++;

			builder.mult = Integer.parseInt(params.get(i).split("=")[1]);
			i++;

			builder.index = moleculei;

			// Atoms
			ArrayList<Atom> atomsL = new ArrayList<>();

			while (!params.get(i).equals("---") && !params.get(i).equals("EXPGEOM")) {
				addAtom(atomsL, params, i);
				i++;
			}
			Atom[] atoms = new Atom[atomsL.size()];
			for (int p = 0; p < atomsL.size(); p++) atoms[p] = atomsL.get(p);
			builder.atoms = atoms;


			// expGeom
			ArrayList<Atom> expGeomL = new ArrayList<>();
			Atom[] expGeom = null;
			if (params.get(i).equals("EXPGEOM")) {
				i++;
				while (!params.get(i).equals("---")) {
					addAtom(expGeomL, params, i);
					i++;
				}
				expGeom = new Atom[expGeomL.size()];
				for (int p = 0; p < expGeomL.size(); p++)
					expGeom[p] = expGeomL.get(p);
			}
			builder.expGeom = expGeom;

			if (expGeomL.size() != 0 && atomsL.size() != expGeomL.size())
				throw new IllegalArgumentException("Atom and expGeom size mismatch!");


			// Datum
			double[] datum = new double[3];

			datum[0] = Double.parseDouble(pncsv.get(datumi).split(" ")[1]);
			datumi++;

			String[] ss = pncsv.get(datumi).split(" ");
			if (ss.length > 1) datum[1] = Double.parseDouble(ss[1]);
			datumi++;

			ss = pncsv.get(datumi).split(" ");
			if (ss.length > 1) datum[2] = Double.parseDouble(ss[1]);
			datumi += 2;

			builder.datum = datum;


			moleculesL.add(builder.build(neededParamsMap));
			i++;

			moleculei++;
		}

		ri.molecules = new RunnableMolecule[moleculesL.size()];
		for (int p = 0; p < moleculesL.size(); p++) {
			ri.molecules[p] = moleculesL.get(p);
		}
		ri.nMolecules = ri.molecules.length;

		outputJSON(ri, "input");
	}
}
