import mndoparam.mndo.MNDOAtom;
import optimize.ParamOptimizer;
import runcycle.AbstractMoleculeRun;
import runcycle.MoleculeRun;
import runcycle.MoleculeRunUnrestricted;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.apache.commons.lang.time.StopWatch;

import datum.DipoleData;
import datum.GeometricalData;
import datum.HeatData;
import datum.IonizationData;

public class Main {
	
  static String trainingset = "CHN";

  public static void main(String[] args) {
	  StopWatch sw = new StopWatch();
	  sw.start();
	  System.err.println ("MNDO Parameterization, updated 27 October. CHN training set (PM7)");
	  
	  boolean useHessian = false;

    for (int numRuns = 0; numRuns < 1; numRuns++) {
    	
    	if (numRuns % 2 == 0) {
    		useHessian = true;
    	}
    	else {
    		useHessian = false;
    	}
    	
    	File input = new File("input.txt");

        File measurements = new File("measurements.txt");
        
        if (useHessian) {
        	measurements = new File("measurementshessian.txt");
        }

        File reference = new File("reference.txt");

        File params = new File("mndoparams.txt");

        try {

          Scanner scanparams = new Scanner(params);

          String[] strings = scanparams.nextLine().split(",");

          double[] paramvector = new double[strings.length];

          for (int i = 0; i < strings.length; i++) {
            paramvector[i] = Double.parseDouble(strings[i]);
          }

          double[] Hparams = new double[13];

          for (int i = 0; i < 13; i++) {
            Hparams[i] = paramvector[i];
          }

          double[] Cparams = new double[13];

          for (int i = 0; i < 13; i++) {
            Cparams[i] = paramvector[i + 13];
          }
          
          double[] Nparams = new double[13];
          
          if (trainingset.contains("N")) {

              for (int i = 0; i < 13; i++) {
                Nparams[i] = paramvector[i + 26];
              }
          }

          File output = new File("mndooutput.txt");

          PrintWriter pw = new PrintWriter(new FileOutputStream(output, true));

          boolean keepgoing = true;

          Scanner sc = new Scanner(input);

          Scanner measurementscan = new Scanner(measurements);

          Scanner referencescan = new Scanner(reference);

          List<ComputationRequest> computationRequests = new ArrayList<>();

          int counter = 0;
          while (keepgoing) {

            ArrayList<MNDOAtom> array = new ArrayList<MNDOAtom>();

            boolean bool;

            int charge;

            int mult;

            if (sc.nextLine().equals("RHF")) {
              bool = true;
              charge = Integer.parseInt(sc.nextLine().substring(7));
              sc.nextLine();
              mult = 1;
            } else {
              bool = false;
              charge = Integer.parseInt(sc.nextLine().substring(7));
              mult = Integer.parseInt(sc.nextLine().substring(5));
            }

            boolean hasGeom = false;

            while (sc.hasNext()) {

              String s = sc.nextLine();
              
              //System.err.println (s);

              if (s.equals("---")) {
                break;
              } else if (s.equals("EXPGEOM")) {
                hasGeom = true;
                break;
              }

              StringTokenizer t = new StringTokenizer(s, " ");
              t.nextToken();
              String element = t.nextToken();
              double d1 = Double.parseDouble(t.nextToken()) * 1.88973;
              double d2 = Double.parseDouble(t.nextToken()) * 1.88973;
              double d3 = Double.parseDouble(t.nextToken()) * 1.88973;

              int Z = 0;

              double[] paramvec = new double[13];

              switch (element) {
                case "H":
                  Z = 1;
                  paramvec = Hparams;
                  break;
                case "C":
                  Z = 6;
                  paramvec = Cparams;
                  break;
                case "N":
                  Z = 7;
                  paramvec = Nparams;
                  break;
                case "O":
                  Z = 8;
                  break;
                case "F":
                  Z = 9;
              }

              MNDOAtom a = new MNDOAtom(new double[]{d1, d2, d3}, Z, paramvec);

              array.add(a);
            }

            MNDOAtom[] atoms = new MNDOAtom[array.size()];

            for (int i = 0; i < array.size(); i++) {
              atoms[i] = array.get(i);

            }


            MNDOAtom[] expgeom = new MNDOAtom[array.size()];

            array.clear();

            if (hasGeom) {
              while (sc.hasNext()) {

                String s = sc.nextLine();

                if (s.equals("---")) {
                  break;
                }


                StringTokenizer t = new StringTokenizer(s, " ");
                t.nextToken();
                String element = t.nextToken();
                double d1 = Double.parseDouble(t.nextToken()) * 1.88973;
                double d2 = Double.parseDouble(t.nextToken()) * 1.88973;
                double d3 = Double.parseDouble(t.nextToken()) * 1.88973;

                int Z = 0;

                double[] paramvec = new double[13];

                switch (element) {
                  case "H":
                    Z = 1;
                    paramvec = Hparams;
                    break;
                  case "C":
                    Z = 6;
                    paramvec = Cparams;
                    break;
                  case "N":
                    Z = 7;
                    paramvec = Nparams;
                    break;
                  case "O":
                    Z = 8;
                    break;
                  case "F":
                    Z = 9;
                }


                MNDOAtom a = new MNDOAtom(new double[]{d1, d2, d3}, Z, paramvec);

                array.add(a);
              }
              for (int i = 0; i < array.size(); i++) {
                expgeom[i] = array.get(i);

              }


            } else { expgeom = null;}

            double[] datum = new double[3];

            boolean hasHessian = false;

            if (measurementscan.nextLine().equals("HEAT HESSIAN")) {
              hasHessian = true;
            }

            datum[0] = Double.parseDouble(referencescan.nextLine().substring(5));

            if (measurementscan.nextLine().equals("DIPOLE")) {
              datum[1] = Double.parseDouble(referencescan.nextLine().substring(7));
            } else {
              referencescan.nextLine();
            }

            if (measurementscan.nextLine().equals("IONIZATION")) {
              datum[2] = Double.parseDouble(referencescan.nextLine().substring(11));
            } else {
              referencescan.nextLine();
            }

            if (referencescan.hasNext()) {
              referencescan.nextLine();
              measurementscan.nextLine();
            }

            computationRequests
                    .add(new ComputationRequest(bool, atoms.clone(), charge, mult, expgeom == null ? null : expgeom
                            .clone(), datum.clone(), hasHessian, counter));
            counter++;
            keepgoing = sc.hasNext();
          }
          
          
          List <ComputationRequest> ParallelComputationRequests = computationRequests.subList(0, 1);
         
          
          int cores = Runtime.getRuntime().availableProcessors();
          pw.println("Running on " + cores + " cores.");
          // Parallel part
          // Run one thread per core (tweakable)
          ForkJoinPool threadPool = new ForkJoinPool(cores);
          // Run parallelMap on thread pool
          List<AbstractMoleculeRun> results = threadPool.submit(() -> ParallelComputationRequests.parallelStream().map(request -> {
            System.out.println("Computation started:" + request.index);
            AbstractMoleculeRun result = request.restricted ? new MoleculeRun(request.atoms, request.charge,
                    request.expgeom, request.datum, request.hasHessian) : new MoleculeRunUnrestricted(request.atoms,
                    request.charge, request.mult, request.expgeom, request.datum, request.hasHessian);

            pw.println("Computation number: " + request.index);
            pw.print("Hessian string: " + result.hessianstr);
            if (result.hessianstr.equals("")) {
              pw.print("\n");
            }
            pw.print("Excel string: " + result.output[0]);
            if (result.output[0].equals("")) {
              pw.print("\n");
            }
            pw.print("HoF string: " + result.output[1]);
            if (result.output[1].equals("")) {
              pw.print("\n");
            }
            pw.print("Dipole string: " + result.output[2]);
            if (result.output[2].equals("")) {
              pw.print("\n");
            }
            pw.print("Ionization string: " + result.output[3]);
            if (result.output[3].equals("")) {
              pw.print("\n");
            }
            pw.print("Geometry string: " + result.output[4]);
            if (result.output[4].equals("")) {
              pw.print("\n");
            }
            pw.print("---------\n");
            // Save to file
            pw.flush();
            System.out.println("Computation complete:" + request.index);
            return result;
          })).get().collect(Collectors.toList());

          // Get output
          List<String[]> outputValues = results.stream().map(result -> result.output.clone()).collect(Collectors.toList());
          List<String> hessianValues = results.stream().map(result -> result.hessianstr).collect(Collectors.toList());
          List<String> optgeoms = results.stream().map(result -> result.newgeomcoords).collect(Collectors.toList());
          
          for (ComputationRequest request : computationRequests.subList(1, computationRequests.size())) {
              System.out.println("Computation started: " + request.index);
              AbstractMoleculeRun result = request.restricted ? new MoleculeRun(request.atoms, request.charge,
                      request.expgeom, request.datum, request.hasHessian) : new MoleculeRunUnrestricted(request.atoms,
                      request.charge, request.mult, request.expgeom, request.datum, request.hasHessian);
              hessianValues.add(result.hessianstr);
              outputValues.add(result.output.clone());
              optgeoms.add (result.newgeomcoords);
              
              pw.println("Computation number: " + request.index);
              pw.print("Hessian string: " + result.hessianstr);
              if (result.hessianstr.equals("")) {
                pw.print("\n");
              }
              pw.print("Excel string: " + result.output[0]);
              if (result.output[0].equals("")) {
                pw.print("\n");
              }
              pw.print("HoF string: " + result.output[1]);
              if (result.output[1].equals("")) {
                pw.print("\n");
              }
              pw.print("Dipole string: " + result.output[2]);
              if (result.output[2].equals("")) {
                pw.print("\n");
              }
              pw.print("Ionization string: " + result.output[3]);
              if (result.output[3].equals("")) {
                pw.print("\n");
              }
              pw.print("Geometry string: " + result.output[4]);
              if (result.output[4].equals("")) {
                pw.print("\n");
              }
              pw.print("---------\n");
              // Save to file
              pw.flush();
              System.out.println("Computation complete:" + request.index);
          }
          
          PrintWriter pw2 = new PrintWriter (new File ("testoutput.txt"));
          
          
          for (int i = 0; i < optgeoms.size(); i++) {
        	  if (i != optgeoms.size() -1) {
        		  pw2.println(optgeoms.get(i) + "---");
        	  }
        	  else {
        		  pw2.print(optgeoms.get(i));
        	  }
          }
          pw2.close();

          pw.print("TOTAL EXCEL STRING:\n");
          for (String[] i : outputValues) {
              pw.print(i[0]);
          }

          pw.print("TOTAL HESSIAN STRING:\n");
          for (String i : hessianValues) {
            pw.print(i);
          }

          pw.print("TOTAL HEAT STRING:\n");
          for (String[] i : outputValues) {
            pw.print(i[1]);
          }

          pw.print("TOTAL DIPOLE STRING:\n");
          for (String[] i : outputValues) {
            pw.print(i[2]);
          }

          pw.print("TOTAL IONIZATION STRING:\n");
          for (String[] i : outputValues) {
            pw.print(i[3]);
          }

          pw.print("TOTAL GEOMETRY STRING:\n");
          for (String[] i : outputValues) {
            pw.print(i[4]);
          }
          
          pw.flush();
          
          
          pw.println("-----------");
          
          int size = 13;
    	  
    	  if (trainingset.equals("CHN")) {
    		  size = 21;
    	  }
          
          double[] sum = new double[size + 1];
          
          for (String[] i: outputValues) {
        	  String[] excelstrsplit = i[0].split(",");
        	  
        	  for (int num = 0; num < size + 1; num++) {
        		  sum[num] += Double.parseDouble(excelstrsplit[num]);
        	  }
          }
          
          
          
          String processedexcelstr = Arrays.toString(Arrays.copyOfRange(sum, 1, size + 1)).substring(1, Arrays.toString(Arrays.copyOfRange(sum, 1, size + 1)).length()-1);
          
          pw.println ("SUM OF ERROR FUNCTION: " + sum[0]);
          pw.println ("GRADIENT SUM: " + processedexcelstr);
          
          pw.flush();
          
          Scanner scan = new Scanner (new File ("MNDOHessianUpdateData.txt"));
          
          String[] strs = null;
    		
          strs = scan.nextLine().substring(15).split(",");
    		
    	  double[] nums = new double[strs.length];
    		
    	  for (int i = 0; i < strs.length; i++) {
    		  nums[i] = Double.parseDouble(strs[i]);
    	  }
    		
    	  int count = 1;
    		
    	  int index = 0;
    	  
    	  
    		
    	  DoubleMatrix hessian = new DoubleMatrix (size, size);
    		
    	  while (index < ((size + 1) * size) / 2) {
    		  for (int i = count - 1; i < size; i++) {
    			  hessian.put(count - 1, i, nums[index]);
    			  hessian.put(i, count - 1, nums[index]);
    			  index++;
    		  }
    		  count++;
    	  }
    		
    		
    	  strs = scan.nextLine().substring(15).split(",");
    	
    	  double[] grad = new double[strs.length];
    		
    	  for (int i = 0; i < strs.length; i++) {
    		  grad[i] = Double.parseDouble(strs[i]);
    	  }
    		
    	  DoubleMatrix oldgradient = new DoubleMatrix (grad);
    		
    		
    	  grad = new double[strs.length];
    		
    	  for (int i = 0; i < strs.length; i++) {
    		  grad[i] = sum[i + 1];
    	  }
    		
    	  DoubleMatrix newgradient = new DoubleMatrix (grad);
    		
    	  strs = scan.nextLine().substring(15).split(",");
    		
    	  double[] dir = new double[strs.length];
    		
    	  for (int i = 0; i < strs.length; i++) {
    		  dir[i] = Double.parseDouble(strs[i]);
    	  }
    		
    	  DoubleMatrix s = new DoubleMatrix (dir);
    	
    	  DoubleMatrix y = newgradient.sub (oldgradient);
    		
    	  double b = y.transpose().mmul(s).get(0);
    		
    		
    	  DoubleMatrix A = y.mmul(y.transpose()).mmul(1/b);
    		
    	  double a = s.transpose().mmul(hessian).mmul(s).get(0);
    		
          DoubleMatrix C = hessian.mmul(s).mmul(s.transpose()).mmul(hessian.transpose()).mmul(1/a);
    		
    	  DoubleMatrix B = hessian.add(A).sub(C);
    	  
    	  if (useHessian) {
    		  double[] hessiansum = new double[(size * (size + 1)) / 2];
    		  for (String hessianstr: hessianValues) {
    			  String[] hessianstrs = hessianstr.split(", ");
    			  for (int i = 0; i < (size * (size + 1)) / 2; i++) {
    				  hessiansum[i] += Double.parseDouble(hessianstrs[i]);
    			  }
    		  }
    		  
    		  index = 0;
    		  count = 1;
    		  while (index < (size * (size + 1)) / 2) {
        		  for (int i = count - 1; i < size; i++) {
        			  B.put(count - 1, i, hessiansum[index]);
        			  B.put(i, count - 1, hessiansum[index]);
        			  index++;
        		  }
        		  count++;
        	  }
    		  
    	  }
    		
    	  String string = "";
    		
          for (int i = 0; i < B.rows; i++) {
    		  for (int j = i; j < B.rows; j++) {
    			  string += B.get(i, j) + ",";
    		  }
    	  }
    	  pw.println ("NEW HESSIAN: " + string);
    	  
    	  pw.println ("HESSIAN EIGENVALUES: " + Eigen.symmetricEigenvalues(B));
    	  
    	  pw.flush();
    	  
    	  pw.close();
    	  
    	 
    	  
    	  ParamOptimizer o = new ParamOptimizer();
    		
    		
    	  for (String[] j: outputValues) {
    		  strs = j[1].split(",");
    			
    		  double[] derivs = new double[strs.length - 2];
    			
    		  for (int i = 2; i < strs.length; i++) {
    			  derivs[i - 2] = Double.parseDouble(strs[i]);
    		  }
    			
    		  o.addData(new HeatData(derivs, Double.parseDouble (strs[0]), Double.parseDouble (strs[1])));
    		  
    		  if (!j[3].equals("")) {
    			  strs = j[3].split(",");
    				
    			  derivs = new double[strs.length - 2];
    				
    			  for (int i = 2; i < strs.length; i++) {
    				  derivs[i - 2] = Double.parseDouble(strs[i]);
    			  }
    				
    			  o.addData(new IonizationData(derivs, Double.parseDouble (strs[0]), Double.parseDouble (strs[1])));
    		  }
    		  
    		  if (!j[2].equals("")) {
    			  strs = j[2].split(",");
    				
    			  derivs = new double[strs.length - 2];
    				
    			  for (int i = 2; i < strs.length; i++) {
    				  derivs[i - 2] = Double.parseDouble(strs[i]);
    			  }
    				
    			  o.addData(new DipoleData(derivs, Double.parseDouble (strs[0]), Double.parseDouble (strs[1])));
    		  }
    		  
    		  if (!j[4].equals("")) {
    			  strs = j[4].split(",");
    				
    			  derivs = new double[strs.length - 2];
    				
    			  for (int i = 2; i < strs.length; i++) {
    				  derivs[i - 2] = Double.parseDouble(strs[i]);
    			  }
    				
    			  o.addData(new GeometricalData(derivs, Double.parseDouble (strs[0]), Double.parseDouble (strs[1])));
    		  }
    		  
    		  
    	  }
    		
    		
    		
    	  o.optimize(B, newgradient);
    	  
    	  PrintWriter write = new PrintWriter(new FileOutputStream(new File("MNDOHessianUpdateData.txt")));
    	  
    	  //write.println("OLD HESSIAN:   " + string);
    	  //write.println("OLD GRADIENT:  " + processedexcelstr);
    	  //write.println("SEARCH VECTOR: " + Arrays.toString(o.changes).substring(1,  Arrays.toString(o.changes).length()-1));
    	  
    	  write.close();
    	  
    	  double[] paramchange = new double[paramvector.length];
          paramchange[0] = o.changes[0];//hardcoded sorry
          paramchange[1] = o.changes[1];
          paramchange[3] = o.changes[2];
          paramchange[5] = o.changes[3];
          paramchange[7] = o.changes[4];
          paramchange[13] = o.changes[5];
          paramchange[14] = o.changes[6];
          paramchange[15] = o.changes[7];
          paramchange[16] = o.changes[8];
          paramchange[17] = o.changes[9];
          paramchange[18] = o.changes[10];
          paramchange[19] = o.changes[11];
          paramchange[20] = o.changes[12];
          
          if (trainingset.contains("N")) {
        	  paramchange[26] = o.changes[13];
              paramchange[27] = o.changes[14];
              paramchange[28] = o.changes[15];
              paramchange[29] = o.changes[16];
              paramchange[30] = o.changes[17];
              paramchange[31] = o.changes[18];
              paramchange[32] = o.changes[19];
              paramchange[33] = o.changes[20];
          }
          for (int num = 0; num < paramvector.length; num++) {
              paramvector[num] = paramvector[num] + paramchange[num];
          }
          
         // write = new PrintWriter(new FileOutputStream(new File("mndoparams.txt")));

          //write.println(Arrays.toString(paramvector).substring(1,  Arrays.toString(paramvector).length()-1));

          //write.close();
          
          write = new PrintWriter(new FileOutputStream(new File("mndooutput.txt"),true));

          write.println(Arrays.toString(paramvector).substring(1,  Arrays.toString(paramvector).length()-1));

          write.close();
    		
          
          
          sw.stop();
          System.out.println("Time " + sw.getTime());


        } catch (Exception e) {
          e.printStackTrace();
        }
    }

  }

}

