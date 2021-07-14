import datum.DipoleData;
import datum.GeometricalData;
import datum.HeatData;
import datum.IonizationData;
import nddoparam.NDDOParams;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import optimize.ParamOptimizer;
import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import runcycle.MoleculeRun;
import runcycle.MoleculeRunRestricted;
import runcycle.MoleculeRunUnrestricted;
import runcycle.input.InputHandler;
import runcycle.input.RawInput;
import runcycle.input.RawMolecule;
import scf.AtomHandler;
import scf.Utils;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

import static scf.Utils.bohr;

public class Main {

    static String trainingSet;
    static int[] atomTypes;
    static ArrayList<String[]> previousOutput = new ArrayList<>(300);

    public static void main(String[] args) {
        StopWatch sw = new StopWatch();
        sw.start();
//        System.out.close();
        AtomHandler.populateAtoms();

        for (int numRuns = 0; numRuns < 1; numRuns++) {
            boolean useHessian = numRuns % 2 == 0; // Hessian every other run
            InputHandler.processInput("input.json");
            RawInput ri = InputHandler.ri;
            System.out.println("MNDO Parameterization, updated 13 July. " + ri.trainingSet + " training set (PM7)");

            NDDOParams[] nddoParams = new MNDOParams[ri.nddoParams.length];
            // TODO change the following line if AM1
            for (int i = 0; i < ri.nddoParams.length; i++) nddoParams[i] = new MNDOParams(ri.nddoParams[i]);

            try {
                List<RawMolecule> requests = new ArrayList<>();
                for (RawMolecule rm : ri.molecules) {
                    requests.add(rm);
                }

                int cores = Runtime.getRuntime().availableProcessors();
                int remainingNonParallel = 5;
                int maxParallel = remainingNonParallel < requests.size() ? requests.size() - remainingNonParallel : 1;
                List<RawMolecule> parallelRequests = requests.subList(0, maxParallel);
                ForkJoinPool threadPool = new ForkJoinPool(cores);

                List<MoleculeRun> results = threadPool.submit(() -> parallelRequests.parallelStream().map(request -> {
                    MoleculeRun result = request.restricted ? new MoleculeRunRestricted(request, nddoParams, ri.atomTypes, useHessian) :
                            new MoleculeRunUnrestricted(request, nddoParams, ri.atomTypes, useHessian);
                    return result;
                })).get().collect(Collectors.toList());

                for (RawMolecule request : requests.subList(maxParallel, requests.size())) {
                    MoleculeRun result = request.restricted ? new MoleculeRunRestricted(request, nddoParams, ri.atomTypes, useHessian) :
                            new MoleculeRunUnrestricted(request, nddoParams, ri.atomTypes, useHessian);
                    results.add(result);
                }

                System.out.println(results.get(0).hessianStr);
//
//                PrintWriter pw2 = new PrintWriter("input.txt");
//                pw2.println("TRAININGSET=" + trainingSet);
//                for (int i = 0; i < outputGeoms.size(); i++) {
//                    if (i != outputGeoms.size() - 1) {
//                        pw2.println(outputGeoms.get(i) + "---");
//                    } else {
//                        pw2.println(outputGeoms.get(i));
//                    }
//                }
//                pw2.close();
//
//                pw.print("TOTAL EXCEL STRING:\n");
//                for (String[] i : outputValues) {
//                    pw.println(i[0]);
//                }
//
//                pw.print("TOTAL HESSIAN STRING:\n");
//                for (String i : hessianValues) {
//                    pw.println(i);
//                }
//
//                pw.print("TOTAL HEAT STRING:\n");
//                for (String[] i : outputValues) {
//                    pw.println(i[1]);
//                }
//
//                pw.print("TOTAL DIPOLE STRING:\n");
//                for (String[] i : outputValues) {
//                    pw.println(i[2]);
//                }
//
//                pw.print("TOTAL IONIZATION STRING:\n");
//                for (String[] i : outputValues) {
//                    pw.println(i[3]);
//                }
//
//                pw.print("TOTAL GEOMETRY STRING:\n");
//                for (String[] i : outputValues) {
//                    pw.println(i[4]);
//                }
//
//                pw.flush();
//
//                pw.println("-----------");
//
//                int size = Utils.getTrainingSetSize(trainingSet);
//                int maxIndex = ((size + 1) * size) / 2;
//
//                double[] sum = new double[size + 1];
//                for (String[] i : outputValues) {
//                    String[] excelstrsplit = i[0].split(",");
//                    for (int num = 0; num < size + 1; num++) {
//                        sum[num] += Double.parseDouble(excelstrsplit[num]);
//                    }
//                }
//                String processedexcelstr = Arrays.toString(Arrays.copyOfRange(sum, 1, size + 1)).substring(1, Arrays.toString(Arrays.copyOfRange(sum, 1, size + 1)).length() - 1);
//
//                pw.println("SUM OF ERROR FUNCTION: " + sum[0]);
//                pw.println("GRADIENT SUM: " + processedexcelstr);
//                pw.flush();
//
//                Scanner scan = new Scanner(new File("MNDOHessianUpdateData.txt"));
//
//                double[] nums = getHessianUpdateData(scan);
//                DoubleMatrix hessian = new DoubleMatrix(size, size);
//                int count = 1;
//                int index = 0;
//                while (index < maxIndex) {
//                    for (int i = count - 1; i < size; i++) {
//                        hessian.put(count - 1, i, nums[index]);
//                        hessian.put(i, count - 1, nums[index]);
//                        index++;
//                    }
//                    count++;
//                }
//
//
//                double[] grad = getHessianUpdateData(scan);
//                DoubleMatrix oldGradient = new DoubleMatrix(grad);
//
//                grad = new double[grad.length];
//                System.arraycopy(sum, 1, grad, 0, grad.length);
//                DoubleMatrix newGradient = new DoubleMatrix(grad);
//
//                double[] dir = getHessianUpdateData(scan);
//                DoubleMatrix s = new DoubleMatrix(dir);
//
//                DoubleMatrix y = newGradient.sub(oldGradient);
//
//                double b = y.transpose().mmul(s).get(0);
//
//
//                DoubleMatrix A = y.mmul(y.transpose()).mmul(1 / b);
//
//                double a = s.transpose().mmul(hessian).mmul(s).get(0);
//
//                DoubleMatrix C = hessian.mmul(s).mmul(s.transpose()).mmul(hessian.transpose()).mmul(1 / a);
//
//                DoubleMatrix B = hessian.add(A).sub(C);
//
//                if (useHessian) {
//                    double[] hessianSum = new double[maxIndex];
//                    for (String hessianStr : hessianValues) {
//                        String[] hessianStrs = hessianStr.split(",");
//                        for (int i = 0; i < maxIndex; i++) {
//                            hessianSum[i] += Double.parseDouble(hessianStrs[i]);
//                        }
//                    }
//
//                    index = 0;
//                    count = 1;
//                    while (index < maxIndex) {
//                        for (int i = count - 1; i < size; i++) {
//                            B.put(count - 1, i, hessianSum[index]);
//                            B.put(i, count - 1, hessianSum[index]);
//                            index++;
//                        }
//                        count++;
//                    }
//                }
//
//                StringBuilder string = new StringBuilder();
//                for (int i = 0; i < B.rows; i++) {
//                    for (int j = i; j < B.rows; j++) {
//                        string.append(B.get(i, j)).append(",");
//                    }
//                }
//                pw.println("NEW HESSIAN: " + string);
//                pw.println("HESSIAN EIGENVALUES: " + Eigen.symmetricEigenvalues(B));
//                pw.flush();
//
//
//                ParamOptimizer o = new ParamOptimizer();
//                for (String[] j : outputValues) { // each j is a molecule
//                    String[] strs = j[1].strip().split(",");
//
//                    double[] derivs = new double[strs.length - 2];
//
//                    for (int i = 2; i < strs.length; i++) {
//                        derivs[i - 2] = Double.parseDouble(strs[i]);
//                    }
//
//                    o.addData(new HeatData(derivs, Double.parseDouble(strs[0]), Double.parseDouble(strs[1])));
//
//                    if (!j[3].equals("")) {
//                        strs = j[3].strip().split(",");
//
//                        derivs = new double[strs.length - 2];
//
//                        for (int i = 2; i < strs.length; i++) {
//                            derivs[i - 2] = Double.parseDouble(strs[i]);
//                        }
//
//                        o.addData(new IonizationData(derivs, Double.parseDouble(strs[0]), Double.parseDouble(strs[1])));
//                    }
//
//                    if (!j[2].equals("")) {
//                        strs = j[2].strip().split(",");
//
//                        derivs = new double[strs.length - 2];
//
//                        for (int i = 2; i < strs.length; i++) {
//                            derivs[i - 2] = Double.parseDouble(strs[i]);
//                        }
//
//                        o.addData(new DipoleData(derivs, Double.parseDouble(strs[0]), Double.parseDouble(strs[1])));
//                    }
//
//                    if (!j[4].equals("")) {
//                        strs = j[4].strip().split(",");
//
//                        derivs = new double[strs.length - 2];
//
//                        for (int i = 2; i < strs.length; i++) {
//                            derivs[i - 2] = Double.parseDouble(strs[i]);
//                        }
//
//                        o.addData(new GeometricalData(derivs, Double.parseDouble(strs[0]), Double.parseDouble(strs[1])));
//                    }
//                }
//                //o.optimize(B, newGradient);
//
//                //PrintWriter write = new PrintWriter(new FileOutputStream(new File("MNDOHessianUpdateData.txt")));
//
//                //write.println("OLD HESSIAN:   " + string);
//                //write.println("OLD GRADIENT:  " + processedexcelstr);
//                //write.println("SEARCH VECTOR: " + Arrays.toString(o.changes).substring(1, Arrays.toString(o.changes).length() - 1));
//
//                //write.close();
//
//                double[] paramChange = new double[paramVector.length];
//                //paramChange[0] = o.changes[0];
//                //paramChange[1] = o.changes[1];
//                //paramChange[3] = o.changes[2];
//                //paramChange[5] = o.changes[3];
//                //paramChange[7] = o.changes[4];
//                String trainingSetNoH = trainingSet.replace("H", "");
//                int startingIndex2 = 5;
//                int z = 8;
//                for (int p = 0; p < trainingSetNoH.length(); p++) {
//                    for (int x = startingIndex2; x < startingIndex2 + 8; x++) {
//                        //paramChange[x + z] = o.changes[x];
//                    }
//                    z += 5;
//                    startingIndex2 += 8;
//                }
//
//                for (int num = 0; num < paramVector.length; num++) {
//                    paramVector[num] = paramVector[num] + paramChange[num];
//                }
//
//                // used to be "write", which is exactly the same as pw
//                //pw.println(Arrays.toString(paramVector).substring(1, Arrays.toString(paramVector).length() - 1));
//                pw.close();

            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        sw.stop();
        System.out.println("Time taken: " + sw.getTime());
    }

    private static double[] getHessianUpdateData(Scanner scan) {
        return Utils.toDoubles(scan.nextLine().split(":")[1].split(","));
    }
}