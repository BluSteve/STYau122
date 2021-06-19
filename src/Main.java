import datum.DipoleData;
import datum.GeometricalData;
import datum.HeatData;
import datum.IonizationData;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import optimize.ParamOptimizer;
import org.apache.commons.lang.time.StopWatch;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import runcycle.MoleculeRun;
import runcycle.MoleculeRunRestricted;
import runcycle.MoleculeRunUnrestricted;
import scf.AtomHandler;
import scf.Utils;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

import static scf.Utils.bohr;

public class Main {

    static String trainingSet;
    static String[] atomTypes;
    static ArrayList<String[]> previousOutput = new ArrayList<>(300);

    public static void main(String[] args) {
        StopWatch sw = new StopWatch();
        sw.start();
        System.err.println("MNDO Parameterization, updated 27 October. CHN training set (PM7)");

        boolean useHessian;

        for (int numRuns = 0; numRuns < 1; numRuns++) { // numRuns < 1 then no Hessian for you
            // one for hessian, one without. this code is verbose.
            useHessian = numRuns % 2 == 1;

            File input = new File("input.txt");
            File reference = new File("reference.txt");
            File params = new File("mndoparams.txt");
            File previous = new File("mndooutput.txt");
            AtomHandler.populateAtoms();

            try {
                // this code is ugly. does mndoparams have a fixed size? TODO a lot of things here
                //readOutput(new Scanner(previous));
                Scanner paramScan = new Scanner(params);
                String[] strings = paramScan.nextLine().split(",");
                double[] paramVector = Utils.toDoubles(strings);

                PrintWriter pw = new PrintWriter(new FileOutputStream("mndooutput.txt", true));

                Scanner sc = new Scanner(input);
                Scanner referenceScan = new Scanner(reference);
                trainingSet = sc.nextLine().split("=")[1];
                atomTypes = trainingSet.split("");

                HashMap<String, MNDOParams> mndoParams = new HashMap<>(5);
                int startingIndex = 0;
                int paramLength = 13;
                int atomNumber = 0;
                while (startingIndex + paramLength <= paramVector.length) {
                    mndoParams.put(atomTypes[atomNumber], new MNDOParams(Arrays.copyOfRange(paramVector, startingIndex, startingIndex + paramLength)));
                    startingIndex += paramLength;
                    atomNumber++;
                }
                System.out.println(mndoParams);


                List<ComputationRequest> requests = new ArrayList<>();
                int counter = 0;
                while (sc.hasNext()) { // for every molecule separated by ---
                    boolean isDoneBefore = false;
                    for (String[] a : previousOutput) {
                        if (counter == Integer.parseInt(a[0])) {
                            isDoneBefore = true;
                            break;
                        }
                    }
                    if (isDoneBefore) continue;

                    ArrayList<MNDOAtom> array = new ArrayList<>();
                    boolean bool = false;
                    int charge, mult;

                    if (sc.nextLine().equals("RHF")) {
                        bool = true;
                        charge = Integer.parseInt(sc.nextLine().split("=")[1]);
                        mult = 1;
                        sc.nextLine();
                    } else {
                        charge = Integer.parseInt(sc.nextLine().split("=")[1]);
                        mult = Integer.parseInt(sc.nextLine().split("=")[1]);
                    }

                    boolean hasGeom = false;

                    while (sc.hasNext()) { // for every line below MULT
                        String s = sc.nextLine();

                        if (s.equals("---")) {
                            break;
                        } else if (s.equals("EXPGEOM")) {
                            hasGeom = true;
                            break;
                        }

                        parseInput(mndoParams, array, s);
                    }

                    MNDOAtom[] atoms = new MNDOAtom[array.size()];
                    atoms = array.toArray(atoms);
                    MNDOAtom[] expGeom = null;
                    array.clear();

                    if (hasGeom) { // for every line below EXPGEOM
                        while (sc.hasNext()) {
                            String s = sc.nextLine();
                            if (s.equals("---")) {
                                break;
                            }
                            parseInput(mndoParams, array, s);
                        }
                        expGeom = new MNDOAtom[array.size()];
                        expGeom = array.toArray(expGeom);
                    }


                    // uses measurements to check if need dipole, ionization, etc. uses values from reference.
                    double[] data = new double[3];
                    data[0] = Double.parseDouble(referenceScan.nextLine().split(" ")[1]); // I wish we standardized delimiters.
                    String[] dipoles = referenceScan.nextLine().split(" ");
                    String kind = "hf_only";
                    if (dipoles.length > 1) {
                        data[1] = Double.parseDouble(dipoles[1]);
                        kind = "limited";
                    }
                    String[] ionizations = referenceScan.nextLine().split(" ");
                    if (ionizations.length > 1) {
                        data[2] = Double.parseDouble(ionizations[1]);
                        kind = "complementary";
                    }
                    if (referenceScan.hasNext()) referenceScan.nextLine();


                    requests
                            .add(new ComputationRequest(bool, atoms.clone(), charge, mult, expGeom == null ? null :
                                    expGeom.clone(), data.clone(), useHessian, kind, counter));
                    counter++;
                }


                // requests contains 1 request for every molecule.
                int cores = Runtime.getRuntime().availableProcessors();
                pw.println("Running on " + cores + " cores.");
                int remainingNonParallel = 14;
                // if requests less than remainingNonParallel then just use parallel computation for one of them
                int maxParallel = remainingNonParallel < requests.size() ? requests.size() - remainingNonParallel : 1;
                List<ComputationRequest> ParallelComputationRequests = requests.subList(0, maxParallel);
                // Parallel part
                // Run one thread per core (tweakable)
                ForkJoinPool threadPool = new ForkJoinPool(cores);
                // Run parallelMap on thread pool
                List<MoleculeRun> results = threadPool.submit(() -> ParallelComputationRequests.parallelStream().map(request -> {
                    MoleculeRun result = request.restricted ? new MoleculeRunRestricted(request.atoms, request.charge,
                            request.expgeom, request.datum, request.hasHessian, trainingSet) : new MoleculeRunUnrestricted(request.atoms,
                            request.charge, request.mult, request.expgeom, request.datum, request.hasHessian, trainingSet);
                    writeOutput(pw, request, result);
                    return result;
                })).get().collect(Collectors.toList());

                // Get output

                List<String[]> outputValues = new ArrayList<>();
                for (String[] slist : previousOutput) {
                    outputValues.add(Arrays.copyOfRange(slist, 2, 6));
                }
                outputValues.addAll(results.stream().map(result -> result.output.clone()).collect(Collectors.toList()));
                // TODO sort this


                List<String> hessianValues = new ArrayList<>();
                for (String[] slist : previousOutput) {
                    hessianValues.add(slist[1]);
                }
                hessianValues.addAll(results.stream().map(result -> result.hessianStr).collect(Collectors.toList()));

                List<String> outputGeoms = new ArrayList<>();
                for (String[] slist : previousOutput) {
                    outputGeoms.add(slist[6]);
                }
                outputGeoms.addAll(results.stream().map(result -> result.newGeomCoords).collect(Collectors.toList()));

                for (ComputationRequest request : requests.subList(maxParallel, requests.size())) {
                    MoleculeRun result = request.restricted ? new MoleculeRunRestricted(request.atoms, request.charge,
                            request.expgeom, request.datum, request.hasHessian, trainingSet) : new MoleculeRunUnrestricted(request.atoms,
                            request.charge, request.mult, request.expgeom, request.datum, request.hasHessian, trainingSet);
                    hessianValues.add(result.hessianStr);
                    outputValues.add(result.output.clone());
                    outputGeoms.add(result.newGeomCoords);

                    writeOutput(pw, request, result);
                }

                PrintWriter pw2 = new PrintWriter("input.txt");
                pw2.println("TRAININGSET=" + trainingSet);
                for (int i = 0; i < outputGeoms.size(); i++) {
                    if (i != outputGeoms.size() - 1) {
                        pw2.println(outputGeoms.get(i) + "---");
                    } else {
                        pw2.println(outputGeoms.get(i));
                    }
                }
                pw2.close();

                pw.print("TOTAL EXCEL STRING:\n");
                for (String[] i : outputValues) {
                    pw.println(i[0]);
                }

                pw.print("TOTAL HESSIAN STRING:\n");
                for (String i : hessianValues) {
                    pw.println(i);
                }

                pw.print("TOTAL HEAT STRING:\n");
                for (String[] i : outputValues) {
                    pw.println(i[1]);
                }

                pw.print("TOTAL DIPOLE STRING:\n");
                for (String[] i : outputValues) {
                    pw.println(i[2]);
                }

                pw.print("TOTAL IONIZATION STRING:\n");
                for (String[] i : outputValues) {
                    pw.println(i[3]);
                }

                pw.print("TOTAL GEOMETRY STRING:\n");
                for (String[] i : outputValues) {
                    pw.println(i[4]);
                }

                pw.flush();

                pw.println("-----------");

                int size = Utils.getTrainingSetSize(trainingSet); // TODO hardcoding, perhaps put in text file reference values for various atoms
                int maxIndex = ((size + 1) * size) / 2;

                double[] sum = new double[size + 1];
                for (String[] i : outputValues) {
                    String[] excelstrsplit = i[0].split(",");
                    for (int num = 0; num < size + 1; num++) {
                        sum[num] += Double.parseDouble(excelstrsplit[num]);
                    }
                }
                String processedexcelstr = Arrays.toString(Arrays.copyOfRange(sum, 1, size + 1)).substring(1, Arrays.toString(Arrays.copyOfRange(sum, 1, size + 1)).length() - 1);

                pw.println("SUM OF ERROR FUNCTION: " + sum[0]);
                pw.println("GRADIENT SUM: " + processedexcelstr);
                pw.flush();

                Scanner scan = new Scanner(new File("MNDOHessianUpdateData.txt"));

                double[] nums = getHessianUpdateData(scan);
                DoubleMatrix hessian = new DoubleMatrix(size, size);
                int count = 1;
                int index = 0;
                while (index < maxIndex) {
                    for (int i = count - 1; i < size; i++) {
                        hessian.put(count - 1, i, nums[index]);
                        hessian.put(i, count - 1, nums[index]);
                        index++;
                    }
                    count++;
                }


                double[] grad = getHessianUpdateData(scan);
                DoubleMatrix oldGradient = new DoubleMatrix(grad);

                grad = new double[grad.length];
                System.arraycopy(sum, 1, grad, 0, grad.length);
                DoubleMatrix newGradient = new DoubleMatrix(grad);

                double[] dir = getHessianUpdateData(scan);
                DoubleMatrix s = new DoubleMatrix(dir);

                DoubleMatrix y = newGradient.sub(oldGradient);

                double b = y.transpose().mmul(s).get(0);


                DoubleMatrix A = y.mmul(y.transpose()).mmul(1 / b);

                double a = s.transpose().mmul(hessian).mmul(s).get(0);

                DoubleMatrix C = hessian.mmul(s).mmul(s.transpose()).mmul(hessian.transpose()).mmul(1 / a);

                DoubleMatrix B = hessian.add(A).sub(C);

                if (useHessian) {
                    double[] hessianSum = new double[maxIndex];
                    for (String hessianStr : hessianValues) {
                        String[] hessianStrs = hessianStr.split(", ");
                        for (int i = 0; i < maxIndex; i++) {
                            hessianSum[i] += Double.parseDouble(hessianStrs[i]);
                        }
                    }

                    index = 0;
                    count = 1;
                    while (index < maxIndex) {
                        for (int i = count - 1; i < size; i++) {
                            B.put(count - 1, i, hessianSum[index]);
                            B.put(i, count - 1, hessianSum[index]);
                            index++;
                        }
                        count++;
                    }
                }

                StringBuilder string = new StringBuilder();
                for (int i = 0; i < B.rows; i++) {
                    for (int j = i; j < B.rows; j++) {
                        string.append(B.get(i, j)).append(",");
                    }
                }
                pw.println("NEW HESSIAN: " + string);
                pw.println("HESSIAN EIGENVALUES: " + Eigen.symmetricEigenvalues(B));
                pw.flush();


                ParamOptimizer o = new ParamOptimizer();
                for (String[] j : outputValues) { // each j is a molecule
                    String[] strs = j[1].strip().split(",");

                    double[] derivs = new double[strs.length - 2];

                    for (int i = 2; i < strs.length; i++) {
                        derivs[i - 2] = Double.parseDouble(strs[i]);
                    }

                    o.addData(new HeatData(derivs, Double.parseDouble(strs[0]), Double.parseDouble(strs[1])));

                    if (!j[3].equals("")) {
                        strs = j[3].strip().split(",");

                        derivs = new double[strs.length - 2];

                        for (int i = 2; i < strs.length; i++) {
                            derivs[i - 2] = Double.parseDouble(strs[i]);
                        }

                        o.addData(new IonizationData(derivs, Double.parseDouble(strs[0]), Double.parseDouble(strs[1])));
                    }

                    if (!j[2].equals("")) {
                        strs = j[2].strip().split(",");

                        derivs = new double[strs.length - 2];

                        for (int i = 2; i < strs.length; i++) {
                            derivs[i - 2] = Double.parseDouble(strs[i]);
                        }

                        o.addData(new DipoleData(derivs, Double.parseDouble(strs[0]), Double.parseDouble(strs[1])));
                    }

                    if (!j[4].equals("")) {
                        strs = j[4].strip().split(",");

                        derivs = new double[strs.length - 2];

                        for (int i = 2; i < strs.length; i++) {
                            derivs[i - 2] = Double.parseDouble(strs[i]);
                        }

                        o.addData(new GeometricalData(derivs, Double.parseDouble(strs[0]), Double.parseDouble(strs[1])));
                    }
                }
                //o.optimize(B, newGradient);

                //PrintWriter write = new PrintWriter(new FileOutputStream(new File("MNDOHessianUpdateData.txt")));

                //write.println("OLD HESSIAN:   " + string);
                //write.println("OLD GRADIENT:  " + processedexcelstr);
                //write.println("SEARCH VECTOR: " + Arrays.toString(o.changes).substring(1, Arrays.toString(o.changes).length() - 1));

                //write.close();

                double[] paramChange = new double[paramVector.length];
                //paramChange[0] = o.changes[0];
                //paramChange[1] = o.changes[1];
                //paramChange[3] = o.changes[2];
                //paramChange[5] = o.changes[3];
                //paramChange[7] = o.changes[4];
                String trainingSetNoH = trainingSet.replace("H", "");
                int startingIndex2 = 5;
                int z = 8;
                for (int p = 0; p < trainingSetNoH.length(); p++) {
                    for (int x = startingIndex2; x < startingIndex2 + 8; x++) {
                        //paramChange[x + z] = o.changes[x];
                    }
                    z += 5;
                    startingIndex2 += 8;
                }

                for (int num = 0; num < paramVector.length; num++) {
                    paramVector[num] = paramVector[num] + paramChange[num];
                }

                // used to be "write", which is exactly the same as pw
                //pw.println(Arrays.toString(paramVector).substring(1, Arrays.toString(paramVector).length() - 1));
                pw.close();

                sw.stop();
                System.out.println("Time taken: " + sw.getTime());
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private static void parseInput(HashMap<String, MNDOParams> mndoParams, ArrayList<MNDOAtom> array, String s) {
        StringTokenizer t = new StringTokenizer(s, " ");
        t.nextToken();
        String element = t.nextToken();
        double x = Double.parseDouble(t.nextToken()) * bohr;
        double y = Double.parseDouble(t.nextToken()) * bohr;
        double z = Double.parseDouble(t.nextToken()) * bohr;

        array.add(new MNDOAtom(AtomHandler.atomsMap.get(element), new double[]{x, y, z}, mndoParams.get(element)));
    }

    private static void writeOutput(PrintWriter pw, ComputationRequest request, MoleculeRun result) {
        System.out.println("Computation started: " + request.index);
        String output = "";
        output += "Computation number: " + request.index + "\n";
        output += "Hessian string: " + result.hessianStr + "\n";
        output += "Excel string: " + result.output[0] + "\n";
        output += "HoF string: " + result.output[1] + "\n";
        output += "Dipole string: " + result.output[2] + "\n";
        output += "Ionization string: " + result.output[3] + "\n";
        output += "Geometry string: " + result.output[4] + "\n";
        output += ("---------\n");

        pw.print(output);
        pw.flush();
        System.out.println("Computation complete: " + request.index);
    }

    private static void readOutput(Scanner s) {
        s.nextLine();
        String[] doubleStrings = new String[7];
        int count = 0;
        while (s.hasNext()) {
            String line = s.nextLine();
            doubleStrings[count] = s.nextLine().strip().split(": ")[1];
            count++;
            if (!line.equals("---------")) {
                previousOutput.add(doubleStrings);
                doubleStrings = new String[7];
                count = 0;
            }
        }
    }

    private static double[] getHessianUpdateData(Scanner scan) {
        return Utils.toDoubles(scan.nextLine().split(":")[1].split(","));
    }
}