package frontend.txt;

import nddo.structs.AtomProperties;
import runcycle.RunIterable;
import runcycle.structs.*;
import tools.Utils;

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
				atomsL.add(toAtom(lines.get(inputi)));
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
			if (ss.length > 1) datum[1] = "DIP=" + ss[1];
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

	private static Atom toAtom(String line) {
		StringTokenizer t = new StringTokenizer(line, " ");
		t.nextToken();
		String name = t.nextToken();
		double[] coords = new double[3];
		for (int q = 0; q < 3; q++) coords[q] = Double.parseDouble(t.nextToken());
		Atom a = new Atom(AtomProperties.getAtomsMap().get(name).getZ(), Utils.bohr(coords));
		return a;
	}

	private static String[] splitCsvLine(String s) {
		return s.split(", *");
	}


	public static void convertFromTXT() throws IOException {
		List<String> pcsv = Files.readAllLines(Path.of("params.csv"));
		int[] Zs = new int[pcsv.size()];
		double[][] params = new double[pcsv.size()][];
		for (int i = 0; i < pcsv.size(); i++) {
			String[] linea = splitCsvLine(pcsv.get(i));
			Zs[i] = AtomProperties.getAtomsMap().get(linea[0]).getZ();
			params[i] = Utils.toDoubles(Arrays.copyOfRange(linea, 1, linea.length));
		}

		List<String> pncsv = Files.readAllLines(Path.of("param-numbers.csv"));
		int[] atomTypes = new int[pncsv.size()];
		int[][] neededParams = new int[pncsv.size()][];
		for (int i = 0; i < pncsv.size(); i++) {
			String[] linea = splitCsvLine(pncsv.get(i));
			atomTypes[i] = AtomProperties.getAtomsMap().get(linea[0]).getZ();
			neededParams[i] = Utils.toInts(Arrays.copyOfRange(linea, 1, linea.length));
		}

		Params p = new Params(Zs, params);
		InputInfo info = new InputInfo(atomTypes, neededParams, p);


		// Molecules
		List<String> mtxt = Files.readAllLines(Path.of("molecules.txt"));

		ArrayList<RunnableMolecule> moleculesL = new ArrayList<>();
		int i = 0;
		int moleculei = 0;

		while (i < mtxt.size()) {
			RunnableMolecule.RMBuilder builder = new RunnableMolecule.RMBuilder();

			String[] minfo = splitCsvLine(mtxt.get(i));
			if (!minfo[0].equals("")) builder.name = minfo[0];

			builder.restricted = minfo[1].equals("RHF");

			builder.charge = Integer.parseInt(minfo[2].split("=")[1]);

			builder.mult = Integer.parseInt(minfo[3].split("=")[1]);

			builder.index = moleculei;

			// Datum
			i++;
			double[] datum = new double[3];
			String[] mdatum = splitCsvLine(mtxt.get(i));

			for (String s : mdatum) {
				if (s.startsWith("HF=")) datum[0] = Double.parseDouble(s.split("=")[1]);
				if (s.startsWith("DIP=")) datum[0] = Double.parseDouble(s.split("=")[1]);
				if (s.startsWith("IE=")) datum[0] = Double.parseDouble(s.split("=")[1]);
			}

			builder.datum = datum;

			// Atoms
			i++;
			ArrayList<Atom> atomsL = new ArrayList<>();

			while (!mtxt.get(i).equals("---") && !mtxt.get(i).equals("EXPGEOM")) {
				atomsL.add(toAtom(mtxt.get(i)));
				i++;
			}
			Atom[] atoms = new Atom[atomsL.size()];
			for (int j = 0; j < atomsL.size(); j++) atoms[j] = atomsL.get(j);
			builder.atoms = atoms;


			// expGeom
			ArrayList<Atom> expGeomL = new ArrayList<>();
			Atom[] expGeom = null;
			if (mtxt.get(i).equals("EXPGEOM")) {
				i++;
				while (!mtxt.get(i).equals("---")) {
					expGeomL.add(toAtom(mtxt.get(i)));
					i++;
				}
				expGeom = new Atom[expGeomL.size()];
				for (int j = 0; j < expGeomL.size(); j++)
					expGeom[j] = expGeomL.get(j);
			}
			builder.expGeom = expGeom;

			if (expGeomL.size() != 0 && atomsL.size() != expGeomL.size())
				throw new IllegalArgumentException("Atom and expGeom size mismatch!");


			moleculesL.add(builder.build(atomTypes, neededParams));
			i++;

			moleculei++;
		}

		RunnableMolecule[] molecules = new RunnableMolecule[moleculesL.size()];
		for (int j = 0; j < moleculesL.size(); j++) {
			molecules[j] = moleculesL.get(j);
		}

		RunIterable iterable = new RunIterable(info, molecules);
		iterable.setLimit(1);
		for (RunOutput ro: iterable) {

		}
	}

	public static void main(String[] args) throws IOException {
		convertFromTXT();
	}
}
