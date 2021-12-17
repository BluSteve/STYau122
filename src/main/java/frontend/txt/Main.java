package frontend.txt;

import nddo.structs.AtomProperties;
import runcycle.IMoleculeResult;
import runcycle.RunIterable;
import runcycle.structs.*;
import tools.Utils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

import static frontend.json.JsonIO.write;

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


			pw.println(String.join(", ", name, rhf, charge, mult));
			pw.println(String.join(", ", datumL));
			pw.write(atoms);
			pw.write(expGeom);
			pw.println("---");


			inputi++;
		}

		pw.close();
	}

	public static RunInput readInput() throws IOException {
		List<String> pcsv = Files.readAllLines(Path.of("params.csv"));
		double[][] params = new double[pcsv.size()][];
		for (int i = 0; i < pcsv.size(); i++) {
			String[] linea = splitCsvLine(pcsv.get(i));
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

		InputInfo info = new InputInfo(atomTypes, neededParams, params);


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
				if (s.startsWith("DIP=")) datum[1] = Double.parseDouble(s.split("=")[1]);
				if (s.startsWith("IE=")) datum[2] = Double.parseDouble(s.split("=")[1]);
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
			if (mtxt.get(i).equals("EXPGEOM")) {
				i++;

				ArrayList<Atom> expGeomL = new ArrayList<>();
				while (!mtxt.get(i).equals("---")) {
					expGeomL.add(toAtom(mtxt.get(i)));
					i++;
				}

				Atom[] expGeom = new Atom[expGeomL.size()];
				for (int j = 0; j < expGeomL.size(); j++)
					expGeom[j] = expGeomL.get(j);

				builder.expGeom = expGeom;
			}


			moleculesL.add(builder.build(atomTypes, neededParams));
			i++;

			moleculei++;
		}

		RunnableMolecule[] molecules = new RunnableMolecule[moleculesL.size()];
		for (int j = 0; j < moleculesL.size(); j++) {
			molecules[j] = moleculesL.get(j);
		}


		return new RunInput(info, molecules);
	}

	public static void updateMolecules(IMoleculeResult[] results) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter("molecules.txt");

		for (IMoleculeResult r : results) {
			RunnableMolecule rm = r.getUpdatedRm();

			pw.println(String.format("%s, %s, CHARGE=%d, MULT=%d", rm.name, rm.restricted ? "RHF" : "UHF",
					rm.charge, rm.mult));

			pw.print(String.format("HF=%.1f", rm.datum[0]));
			if (rm.datum[1] != 0) pw.print(String.format(", DIP=%.1f", rm.datum[1]));
			if (rm.datum[2] != 0) pw.print(String.format(", IE=%.1f", rm.datum[2]));
			pw.println();

			String precision = "%19.16f";
			String format = String.format("%%d    %%s    %s    %s    %s", precision, precision, precision);
			for (int i = 0; i < rm.atoms.length; i++) {
				Atom atom = rm.atoms[i];

				double[] coords = Utils.debohr(atom.coords);

				pw.println(String.format(format, i + 1,
						AtomProperties.getAtoms()[atom.Z].getName(), coords[0], coords[1], coords[2]));
			}

			if (r.isExpAvail()) {
				pw.println("EXPGEOM");

				for (int i = 0; i < rm.expGeom.length; i++) {
					Atom atom = rm.expGeom[i];

					double[] coords = Utils.debohr(atom.coords);

					pw.println(String.format(format, i + 1,
							AtomProperties.getAtoms()[atom.Z].getName(), coords[0], coords[1], coords[2]));
				}
			}

			pw.println("---");
		}

		pw.close();
	}

	public static void updateParams(InputInfo info) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter("params.csv");

		for (int i = 0; i < info.atomTypes.length; i++) {
			String name = AtomProperties.getAtoms()[info.atomTypes[i]].getName();
			String paramsStr = Arrays.toString(info.getParams()[i]);

			paramsStr = paramsStr.substring(1, paramsStr.length() - 1);

			pw.println(String.format("%s, %s", name, paramsStr));
		}

		pw.close();
	}

	public static void updateOutput(RunOutput ro) throws FileNotFoundException {
		updateMolecules(ro.results);
		updateParams(ro.nextRunInfo);
	}

	private static Atom toAtom(String line) {
		StringTokenizer t = new StringTokenizer(line, " ");
		t.nextToken();
		String name = t.nextToken();
		double[] coords = new double[3];
		for (int q = 0; q < 3; q++) coords[q] = Double.parseDouble(t.nextToken());
		return new Atom(AtomProperties.getAtomsMap().get(name).getZ(), Utils.bohr(coords));
	}

	private static String[] splitCsvLine(String s) {
		return s.split(", *");
	}

	public static void main(String[] args) throws IOException {
		Files.createDirectories(Path.of("pastinputs"));
		Files.createDirectories(Path.of("outputs"));

		RunInput input = readInput();

		RunIterable iterable = new RunIterable(input, 10);
		int i = 0;
		for (RunOutput ro : iterable) {
			write(ro.getInput(), String.format("pastinputs/%04d-%s", i, ro.getInput().hash));
			write(ro, String.format("outputs/%04d-%s-%s", i, ro.getInput().hash, ro.hash));
			i++;
		}

		System.exit(0);
	}
}
