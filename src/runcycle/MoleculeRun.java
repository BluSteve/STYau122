package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimization;
import nddoparam.NDDOSolution;
import nddoparam.param.ParamGradientAnalytical;
import nddoparam.param.ParamGradientRestricted;
import nddoparam.param.ParamHessian;
import org.jblas.DoubleMatrix;

public abstract class MoleculeRun {
    protected static double LAMBDA = 1E-7;
    public String[] output;
    public int[] atomTypes;
    public String hessianStr, newGeomCoords;
    protected NDDOAtom[] atoms, expGeom;
    protected int charge, size, mult;
    protected boolean runHessian;
    protected double[] datum;
    protected ParamGradientAnalytical g;
    protected ParamHessian h;
    protected NDDOGeometryOptimization opt;
    protected NDDOSolution expSolution;
    protected String totalHeatDeriv = "";
    protected String totalIonizationDeriv = "";
    protected String totalDipoleDeriv = "";
    protected String totalGeomDeriv = "";
    protected String totalExcelStr = "";
    protected String excelStr = "";
    protected String excelStr2 = "";
    protected String kind;


    public MoleculeRun(NDDOAtom[] atoms, int charge, NDDOAtom[] expGeom, double[] datum, boolean runHessian, String kind, int[] atomTypes, int mult) {
        this.atomTypes = atomTypes;
        this.atoms = atoms;
        this.expGeom = expGeom;
        this.charge = charge;
        this.runHessian = runHessian;
        this.datum = datum; // reference heat + dipole + ionization. size = 1, 2 or 3
        this.size = atomTypes.length;
        this.mult = mult;
        this.kind = kind;
    }


    protected void generateGeomCoords() {
        StringBuilder sb = new StringBuilder(newGeomCoords);
        for (int i = 0; i < atoms.length; i++) {
            String s = (i + 1) + "    " + atoms[i].getName() + "    " + String.format("%.9f", atoms[i].getCoordinates()[0] / 1.88973) + "    " + String.format("%.9f", atoms[i].getCoordinates()[1] / 1.88973) + "    " + String.format("%.9f", atoms[i].getCoordinates()[2] / 1.88973) + "\n";
            sb.append(s);
        }

        if (this.expGeom != null) {
            sb.append("EXPGEOM\n");
            for (int i = 0; i < this.expGeom.length; i++) {
                String s = (i + 1) + "    " + this.expGeom[i].getName() + "    " + String.format("%.9f", this.expGeom[i].getCoordinates()[0] / 1.88973) + "    " + String.format("%.9f", this.expGeom[i].getCoordinates()[1] / 1.88973) + "    " + String.format("%.9f", this.expGeom[i].getCoordinates()[2] / 1.88973) + "\n";
                sb.append(s);
            }
        }
        newGeomCoords = sb.toString();
    }

    protected void runHessian() {
        DoubleMatrix paramIndexes = new DoubleMatrix(size, 2);

        for (int i = 0; i < size; i++) {
            if (i < 5) {
                paramIndexes.put(i, 0, 1);
            } else if (i < 13) {
                paramIndexes.put(i, 0, 6);
            } else if (i < 21) {
                paramIndexes.put(i, 0, 7);
            } else if (i < 29) {
                paramIndexes.put(i, 0, 8);
            } else if (i < 37) {
                paramIndexes.put(i, 0, 9);
            }
        }

        paramIndexes.put(0, 1, 0);
        paramIndexes.put(1, 1, 1);
        paramIndexes.put(2, 1, 3);
        paramIndexes.put(3, 1, 5);
        paramIndexes.put(4, 1, 7);

        for (int i = 0; i < 8; i++) {
            paramIndexes.put(i + 5, 1, i);
        }
        if (size > 13) {
            for (int i = 0; i < 8; i++) {
                paramIndexes.put(i + 13, 1, i);
            }
        }

        StringBuilder hessianSB = new StringBuilder();
        for (int I = 0; I < size; I++) {
            for (int j = I; j < size; j++) {
                // TODO optimize this
                int Z1 = (int) Math.round(paramIndexes.get(I, 0));
                int Z2 = (int) Math.round(paramIndexes.get(j, 0));
                int param1 = (int) Math.round(paramIndexes.get(I, 1));
                int param2 = (int) Math.round(paramIndexes.get(j, 1));

                getH(Z1, param1, Z2, param2);

                h.constructErrors(datum[0]);

                if (this.expGeom != null) {
                    h.createExpGeom(this.expGeom, expSolution);
                    h.addGeomError();
                }


                if (datum[1] != 0) {
                    h.addDipoleError(datum[1]);
                }

                if (datum[2] != 0) {
                    h.addIEError(datum[2]);
                }

                hessianSB.append(h.hessian()).append(",");
            }
        }
        hessianStr = hessianSB.toString();
    }

    protected void runGradient() {
        g = new ParamGradientRestricted(opt.s, kind, datum, expSolution); // time intensive step
        StringBuilder HFDerivsSB = new StringBuilder(datum[0] + "," + g.getS().hf + ",");
        StringBuilder dipoleDerivsSB = new StringBuilder(datum[1] + "," + g.getS().dipole + ",");
        StringBuilder IEDerivsSB = new StringBuilder(datum[2] + "," + -g.getS().homo + ",");
        StringBuilder geomDerivsSB = new StringBuilder("0," + g.getE().geomGradient + ",");
        StringBuilder mainDataSB = new StringBuilder();
        mainDataSB.append(g.getE().getTotalError());

        for (int atomType : this.atomTypes) {
            switch (kind) {
                case "a":
                    appendToSB(g.getHFDerivs()[atomType], HFDerivsSB);
                    break;
                case "b":
                    appendToSB(g.getHFDerivs()[atomType], HFDerivsSB);
                    appendToSB(g.getDipoleDerivs()[atomType], dipoleDerivsSB);
                    break;
                case "c":
                    appendToSB(g.getHFDerivs()[atomType], HFDerivsSB);
                    appendToSB(g.getIEDerivs()[atomType], IEDerivsSB);
                    break;
                case "d":
                    appendToSB(g.getHFDerivs()[atomType], HFDerivsSB);
                    appendToSB(g.getDipoleDerivs()[atomType], dipoleDerivsSB);
                    appendToSB(g.getIEDerivs()[atomType], IEDerivsSB);
                    break;
            }
            if (expSolution != null) appendToSB(g.getGeomDerivs()[atomType], geomDerivsSB);
            appendToSB(g.getTotalDerivs()[atomType], mainDataSB);
        }

        appendToSB(new double[] {datum[0], g.getS().hf, datum[1], g.getS().dipole,
                        datum[2], g.getS().homo, 0, g.getE().geomGradient}, mainDataSB);

        totalHeatDeriv = HFDerivsSB.toString();
        totalDipoleDeriv = dipoleDerivsSB.toString();
        totalIonizationDeriv = IEDerivsSB.toString();
        totalGeomDeriv = geomDerivsSB.toString();
        totalExcelStr = mainDataSB.toString();
    }

    private static void appendToSB(double[] array, StringBuilder sb) {
        for (int i = 0; i < array.length - 1; i++) sb.append(array[i]).append(',');
        sb.append(array[array.length - 1]);
    }

//    protected void runGradient() {
//        StringBuilder totalHeatDerivSB = new StringBuilder();
//        StringBuilder totalDipoleDerivSB = new StringBuilder();
//        StringBuilder totalIonizationDerivSB = new StringBuilder();
//        StringBuilder totalGeomDerivSB = new StringBuilder();
//        StringBuilder excelSB = new StringBuilder();
//        StringBuilder excelSB2 = new StringBuilder();
//        for (int numit = 0; numit < 8; numit++) {
//            if (numit != 2 && numit != 4 && numit != 6) {
//                getG(1, numit, mult); // Z=1, numit is the paramNum. so in this case it would be 0,1,3,5,7
//                if (numit == 0) {
//                    totalHeatDerivSB.append(datum[0]).append(", ").append(g.getEHf()).append(", ");
//                }
//                g.constructErrors(datum[0]);
//                if (this.expGeom != null) {
//                    g.createExpGeom(this.expGeom, expSolution);
//                    g.addGeomError();
//                }
//                if (datum[1] != 0) {
//                    if (numit == 0) { // wait why is there no need to compute the finite difference gradient here
//                        totalDipoleDerivSB.append(datum[1]).append(", ").append(g.getEDipole()).append(", ");
//                    }
//
//                    g.addDipoleError(datum[1]);
//                    totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole())).append(", ");
//                }
//
//                if (datum[2] != 0) {
//                    if (numit == 0) {
//                        totalIonizationDerivSB.append(datum[2]).append(", ").append(-g.getEHomo()).append(", ");
//                    }
//                    g.addIEError(datum[1]);
//                    totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo())).append(", ");
//                }
//
//                if (this.expGeom != null) {
//                    if (numit == 0) {
//                        totalGeomDerivSB.append("0, ").append(g.getEGradient()).append(",");
//                    }
//                    totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient())).append(", ");
//                }
//
//                System.out.println(g.gradient());
//                totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf())).append(", ");
//
//                excelSB.append(",").append(g.gradient());
//            }
//        }
//
//        NDDOSolutionRestricted soln = (NDDOSolutionRestricted) opt.s;
//        g = new ParamGradientRestricted(soln, this.kind);
//
//
//        for (int numit = 0; numit < 8; numit++) {
//            getG(6, numit, mult);
//
//            excelStr2 = "," + datum[0] + "," + g.getEHf() + ",";
//
//
//            g.constructErrors(datum[0]);
//
//            if (this.expGeom != null) {
//                g.createExpGeom(this.expGeom, expSolution);
//                g.addGeomError();
//            }
//
//            if (datum[1] != 0) {
//                g.addDipoleError(datum[1]);
//
//                excelSB2.append(",").append(datum[1]).append(",").append(g.getEDipole()).append(",");
//
//                if (numit == 7 && trainingSet.equals("HC")) {
//                    totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole()));
//                } else {
//                    totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole())).append(", ");
//                }
//
//            } else {
//                excelSB2.append(",,,");
//            }
//
//            if (datum[2] != 0) {
//                g.addIEError(datum[2]);
//
//                excelSB2.append(",").append(datum[2]).append(",").append(-g.getEHomo()).append(",");
//                if (numit == 7 && trainingSet.equals("HC")) {
//                    totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo()));
//                } else {
//                    totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo())).append(", ");
//                }
//            } else {
//                excelSB2.append(",,,");
//            }
//
//            if (this.expGeom != null) {
//                excelSB2.append("," + 0 + ",").append(g.getEGradient()).append(",");
//                if (numit == 7 && trainingSet.equals("HC")) {
//                    totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient()));
//                } else {
//                    totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient())).append(", ");
//                }
//            }
//
//            if (numit == 7 && trainingSet.equals("HC")) {
//                totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf()));
//            } else {
//                totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf())).append(", ");
//            }
//
//
//            if (this.mult == 1 && numit > 0 && numit < 7) {
//                double Hfderiv = 1 / LAMBDA * (g.getEPrimeHf() - g.getEHf());
//
//                double test = Hfderivs[numit - 1];
//
//                if (Math.abs(Hfderiv - test) > 1E-2) {
//                    System.err.println("Something broke. Again. " + numit);
//                    System.err.println("yay");
//                    System.err.println(Hfderiv);
//                    System.err.println(test);
//                    System.exit(0);
//                } else {
//                    System.err.println("the Hf derivatives work");
//                    System.err.println(Hfderiv);
//                    System.err.println(test);
//                }
//
//
//            }
//
//            if (this.mult == 1 && numit > 0 && numit < 7) {
//                double IEderiv = 1 / LAMBDA * (g.getEPrimeHomo() - g.getEHomo());
//
//                double test = IEderivs[numit - 1];
//
//                if (Math.abs(IEderiv - test) > 1E-2) {
//                    System.err.println("Something broke. Again. " + numit);
//                    System.err.println("yay");
//                    System.err.println(IEderiv);
//                    System.err.println(test);
//                    System.exit(0);
//                } else {
//                    System.err.println("the IE derivatives work");
//                    System.err.println(IEderiv);
//                    System.err.println(test);
//                }
//
//
//            }
//
//
//            if (this.mult == 1 && numit > 0 && numit < 7) {
//
//                NDDOSolutionRestricted solp = (NDDOSolutionRestricted) g.getEPrime().soln;
//                NDDOSolutionRestricted sol = (NDDOSolutionRestricted) g.getE().soln;
//
//
//                DoubleMatrix densityderiv = (solp.densityMatrix().sub(sol.densityMatrix()).mmul(1 / LAMBDA));
//
//
//                if (!NDDOSolution.isSimilar(densityderivs[numit - 1], densityderiv, 1E-5)) {
//                    System.err.println("Something broke. Again. " + numit);
//                    System.err.println("yay");
//                    System.err.println(densityderiv);
//                    System.err.println(densityderivs[numit - 1]);
//                    System.exit(0);
//                } else {
//                    System.err.println("the density derivatives work");
//                }
//
//
//            }
//
//            if (this.mult == 1 && numit > 0 && numit < 7) {
//                double Dipolederiv = 1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole());
//
//                double test = Dipolederivs[numit - 1];
//
//                if (Math.abs(Dipolederiv - test) > 1E-2) {
//                    System.err.println("Something broke. Again. " + numit);
//                    System.err.println("yay");
//                    System.err.println(Dipolederiv);
//                    System.err.println(test);
//                    System.exit(0);
//                } else {
//                    System.err.println("the dipole derivatives work");
//                    System.err.println(Dipolederiv);
//                    System.err.println(test);
//                }
//
//
//            }
//
//
//            excelSB.append(",").append(g.gradient());
//
//        }
//
//        if (trainingSet.strip().contains("N")) {
//            for (int numit = 0; numit < 8; numit++) {
//
//                getG(7, numit, mult);
//
//
//                excelStr2 = "," + datum[0] + "," + g.getEHf() + ",";
//
//
//                g.constructErrors(datum[0]);
//
//                if (this.expGeom != null) {
//                    g.createExpGeom(this.expGeom, expSolution);
//                    g.addGeomError();
//                }
//
//                if (datum[1] != 0) {
//                    g.addDipoleError(datum[1]);
//
//                    excelSB2.append(",").append(datum[1]).append(",").append(g.getEDipole()).append(",");
//
//                    if (numit == 7) {
//                        totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole()));
//                    } else {
//                        totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole())).append(", ");
//                    }
//
//                } else {
//                    excelSB2.append(",,,");
//                }
//
//                if (datum[2] != 0) {
//                    g.addIEError(datum[2]);
//
//                    excelSB2.append(",").append(datum[2]).append(",").append(-g.getEHomo()).append(",");
//                    if (numit == 7) {
//                        totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo()));
//                    } else {
//                        totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo())).append(", ");
//                    }
//                } else {
//                    excelSB2.append(",,,");
//                }
//
//                if (this.expGeom != null) {
//                    excelSB2.append("," + 0 + ",").append(g.getEGradient()).append(",");
//                    if (numit == 7) {
//                        totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient()));
//                    } else {
//                        totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient())).append(", ");
//                    }
//                }
//
//
//                if (numit == 7) {
//                    totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf()));
//                } else {
//                    totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf())).append(", ");
//                }
//
//                excelSB.append(",").append(g.gradient());
//            }
//        }
//        totalDipoleDeriv = totalDipoleDerivSB.toString();
//        totalGeomDeriv = totalGeomDerivSB.toString();
//        totalHeatDeriv = totalHeatDerivSB.toString();
//        totalIonizationDeriv = totalIonizationDerivSB.toString();
//        excelStr = excelSB.toString();
//        excelStr2 = excelSB2.toString();
//    }

    protected void outputErrorFunction() {
        System.out.println("Error function: " + g.getE().getTotalError());

//        totalExcelStr += g.getE().getTotalError() + excelStr + excelStr2;

        output = new String[]{totalExcelStr, totalHeatDeriv, totalDipoleDeriv, totalIonizationDeriv, totalGeomDeriv};
    }

    protected abstract void getH(int Z1, int param1, int Z2, int param2);


}