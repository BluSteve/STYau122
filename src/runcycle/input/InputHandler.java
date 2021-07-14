package runcycle.input;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import nddoparam.NDDOAtom;
import scf.AtomHandler;
import scf.AtomProperties;
import scf.Utils;

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
            ri = (new Gson()).fromJson(new FileReader(path), RawInput.class);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void convertFromTXT() {
        AtomHandler.populateAtoms();
        ri = new RawInput();
        try {
            List<String> lines = Files.readAllLines(Path.of("input.txt"));
            List<String> datums = Files.readAllLines(Path.of("reference.txt"));
            ri.trainingSet = lines.get(0).split("=")[1];
            ArrayList<RawMolecule> moleculesL = new ArrayList<>();

            int i = 1;
            try {
                while (i < lines.size()) {
                    RawMolecule rm = new RawMolecule();
                    ArrayList<RawAtom> atomsL = new ArrayList<>();
                    ArrayList<RawAtom> expGeomL = new ArrayList<>();

                    rm.uhf = lines.get(i).equals("UHF");
                    i++;

                    rm.charge = Integer.parseInt(lines.get(i).split("=")[1]);
                    i++;

                    rm.mult = Integer.parseInt(lines.get(i).split("=")[1]);
                    i++;

                    while (!lines.get(i).equals("---") && !lines.get(i).equals("EXPGEOM")) {
                        addAtom(atomsL, lines, i);
                        i++;
                    }

                    StringBuilder nameBuilder = new StringBuilder();
                    HashMap<String, Integer> nameOccurrences = new HashMap<>();
                    ArrayList<Integer> tempZs = new ArrayList<>(Utils.maxAtomNum);
                    for (RawAtom a : atomsL) {
                        AtomProperties ap = AtomHandler.atomsMap.get(a.name);
                        rm.nElectrons += ap.getQ();
                        if (!tempZs.contains(ap.getZ())) {
                            tempZs.add(ap.getZ());
                        }
                        if (!nameOccurrences.containsKey(a.name)) nameOccurrences.put(a.name, 1);
                        else nameOccurrences.put(a.name, nameOccurrences.get(a.name) + 1);

                    }
                    for (String key : nameOccurrences.keySet()) {
                        nameBuilder.append(key).append(nameOccurrences.get(key));
                    }
                    Collections.sort(tempZs);
                    int[] uniqueZs = new int[tempZs.size()];
                    for (int u = 0; u < tempZs.size(); u++) uniqueZs[u] = tempZs.get(u);
                    rm.uniqueZs = uniqueZs;
                    rm.name = nameBuilder.toString();


                    if (lines.get(i).equals("EXPGEOM")) {
                        i++;
                        while (!lines.get(i).equals("---")) {
                            addAtom(expGeomL, lines, i);
                            i++;
                        }
                        rm.expGeom = new RawAtom[atomsL.size()];
                        for (int p = 0; p < expGeomL.size(); p++) rm.expGeom[p] = expGeomL.get(p);
                    } else {
                        rm.expGeom = new RawAtom[0];
                    }

                    rm.atoms = new RawAtom[atomsL.size()];
                    for (int p = 0; p < atomsL.size(); p++) rm.atoms[p] = atomsL.get(p);
                    moleculesL.add(rm);
                    i++;
                }
            } catch (IndexOutOfBoundsException e) { }
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
            }
            ri.molecules = new RawMolecule[moleculesL.size()];
            for (int p = 0; p < moleculesL.size(); p++) ri.molecules[p] = moleculesL.get(p);
            for (int j = 0; j < ri.molecules.length; j++) ri.molecules[j].index = j;

            GsonBuilder builder = new GsonBuilder();
            builder.setPrettyPrinting();
            Gson gson = builder.create();
            FileWriter fw = new FileWriter("input.json");
            gson.toJson(ri, fw);
            fw.close();
        } catch (
                IOException e) {
            e.printStackTrace();
        }

    }

    private static void addAtom(List<RawAtom> rawAtoms, List<String> lines, int i) {
        RawAtom ra = new RawAtom();
        StringTokenizer t = new StringTokenizer(lines.get(i), " ");
        t.nextToken();
        ra.name = t.nextToken();
        for (int q = 0; q < 3; q++) ra.coords[q] = Double.parseDouble(t.nextToken());
        rawAtoms.add(ra);
    }

    public static void main(String[] args) {
        InputHandler.convertFromTXT();
    }
}
