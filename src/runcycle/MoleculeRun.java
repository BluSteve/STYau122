package runcycle;

import nddoparam.NDDOAtom;
import nddoparam.NDDOGeometryOptimization;
import nddoparam.NDDOSolution;
import nddoparam.param.NDDOParamGradient;
import nddoparam.param.NDDOParamHessian;
import org.jblas.DoubleMatrix;
import scf.Utils;

public abstract class MoleculeRun {
    public String[] output;
    public String hessianStr, newGeomCoords, trainingSet;
    protected NDDOAtom[] atoms, expGeom;
    protected int charge, size, mult;
    protected boolean runHessian;
    protected double[] datum;
    protected NDDOParamGradient g;
    protected NDDOParamHessian h;
    protected NDDOGeometryOptimization opt;
    protected NDDOSolution expSolution;
    protected static double LAMBDA = 1E-7;
    protected String totalHeatDeriv = "";
    protected String totalIonizationDeriv = "";
    protected String totalDipoleDeriv = "";
    protected String totalGeomDeriv = "";
    protected String totalExcelStr = "";
    protected String excelStr = "";
    protected String excelStr2 = "";


    public MoleculeRun(NDDOAtom[] atoms, int charge, NDDOAtom[] expGeom, double[] datum, boolean runHessian, String trainingSet, int mult) {
        this.trainingSet = trainingSet;
        this.atoms = atoms;
        this.expGeom = expGeom;
        this.charge = charge;
        this.runHessian = runHessian;
        this.datum = datum;
        this.size = Utils.getTrainingSetSize(trainingSet);
        this.mult = mult;
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

                hessianSB.append(h.hessian()).append(", ");

            }
        }
        hessianStr = hessianSB.toString();
        hessianStr += "\n";
    }

    protected void runGradient() {
        StringBuilder totalHeatDerivSB = new StringBuilder();
        StringBuilder totalDipoleDerivSB =new StringBuilder();
        StringBuilder totalIonizationDerivSB=new StringBuilder();
        StringBuilder totalGeomDerivSB=new StringBuilder();
        StringBuilder excelSB = new StringBuilder();
        StringBuilder excelSB2 = new StringBuilder();
        for (int numit = 0; numit < 8; numit++) {
            if (numit != 2 && numit != 4 && numit != 6) {
                getG(1, numit, mult);

                if (numit == 0) {
                    totalHeatDerivSB.append(datum[0]).append(", ").append(g.getEHf()).append(", ");
                }


                g.constructErrors(datum[0]);

                if (this.expGeom != null) {
                    g.createExpGeom(this.expGeom, expSolution);
                    g.addGeomError();
                }


                if (datum[1] != 0) {
                    if (numit == 0) {
                        totalDipoleDerivSB.append(datum[1]).append(", ").append(g.getEDipole()).append(", ");
                    }

                    g.addDipoleError(datum[1]);
                    totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole())).append(", ");
                }

                if (datum[2] != 0) {
                    if (numit == 0) {
                        totalIonizationDerivSB.append(datum[2]).append(", ").append(-g.getEHomo()).append(", ");
                    }
                    g.addIEError(datum[1]);
                    totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo())).append(", ");
                }

                if (this.expGeom != null) {
                    if (numit == 0) {
                        totalGeomDerivSB.append("0, ").append(g.getEGradient()).append(",");
                    }
                    totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient())).append(", ");
                }
                
                System.out.println(g.gradient());
                totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf())).append(", ");
                excelSB.append(",").append(g.gradient());
            }
        }
        for (int numit = 0; numit < 8; numit++) {
            getG(6, numit, mult);

            excelStr2 = "," + datum[0] + "," + g.getEHf() + ",";


            g.constructErrors(datum[0]);

            if (this.expGeom != null) {
                g.createExpGeom(this.expGeom, expSolution);
                g.addGeomError();
            }

            if (datum[1] != 0) {
                g.addDipoleError(datum[1]);

                excelSB2.append(",").append(datum[1]).append(",").append(g.getEDipole()).append(",");

                if (numit == 7 && trainingSet.equals("CH")) {
                    totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole())).append("\n");
                } else {
                    totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole())).append(", ");
                }

            } else {
                excelSB2.append(",,,");
            }

            if (datum[2] != 0) {
                g.addIEError(datum[2]);

                excelSB2.append(",").append(datum[2]).append(",").append(-g.getEHomo()).append(",");
                if (numit == 7 && trainingSet.equals("CH")) {
                    totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo())).append("\n");
                } else {
                    totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo())).append(", ");
                }
            } else {
                excelSB2.append(",,,");
            }

            if (this.expGeom != null) {
                excelSB2.append("," + 0 + ",").append(g.getEGradient()).append(",");
                if (numit == 7 && trainingSet.equals("CH")) {
                    totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient())).append("\n");
                } else {
                    totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient())).append(", ");
                }
            }

            if (numit == 7 && trainingSet.equals("CH")) {
                totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf())).append("\n");
            } else {
                totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf())).append(", ");
            }

            excelSB.append(",").append(g.gradient());

        }

        if (trainingSet.strip().contains("N")) {
            for (int numit = 0; numit < 8; numit++) {

                getG(7, numit, mult);


                excelStr2 = "," + datum[0] + "," + g.getEHf() + ",";


                g.constructErrors(datum[0]);

                if (this.expGeom != null) {
                    g.createExpGeom(this.expGeom, expSolution);
                    g.addGeomError();
                }

                if (datum[1] != 0) {
                    g.addDipoleError(datum[1]);

                    excelSB2.append(",").append(datum[1]).append(",").append(g.getEDipole()).append(",");

                    if (numit == 7) {
                        totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole())).append("\n");
                    } else {
                        totalDipoleDerivSB.append(1 / LAMBDA * (g.getEPrimeDipole() - g.getEDipole())).append(", ");
                    }

                } else {
                    excelSB2.append(",,,");
                }

                if (datum[2] != 0) {
                    g.addIEError(datum[2]);

                    excelSB2.append(",").append(datum[2]).append(",").append(-g.getEHomo()).append(",");
                    if (numit == 7) {
                        totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo())).append("\n");
                    } else {
                        totalIonizationDerivSB.append(1 / LAMBDA * (-g.getEPrimeHomo() + g.getEHomo())).append(", ");
                    }
                } else {
                    excelSB2.append(",,,");
                }

                if (this.expGeom != null) {
                    excelSB2.append("," + 0 + ",").append(g.getEGradient()).append(",");
                    if (numit == 7) {
                        totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient())).append("\n");
                    } else {
                        totalGeomDerivSB.append(1 / LAMBDA * (g.getEPrimeGradient() - g.getEGradient())).append(", ");
                    }
                }


                if (numit == 7) {
                    totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf())).append("\n");
                } else {
                    totalHeatDerivSB.append(1 / LAMBDA * (g.getEPrimeHf() - g.getEHf())).append(", ");
                }

                excelSB.append(",").append(g.gradient());
            }
        }
        totalDipoleDeriv = totalDipoleDerivSB.toString();
        totalGeomDeriv = totalGeomDerivSB.toString();
        totalHeatDeriv = totalHeatDerivSB.toString();
        totalIonizationDeriv = totalIonizationDerivSB.toString();
        excelStr = excelSB.toString();
        excelStr2 = excelSB2.toString();
    }

    protected void outputErrorFunction() {
        System.out.println("Error function: " + g.getE().constructErrorFunction());

        totalExcelStr += g.getE().constructErrorFunction() + excelStr + excelStr2 + "\n";

        output = new String[]{totalExcelStr, totalHeatDeriv, totalDipoleDeriv, totalIonizationDeriv, totalGeomDeriv};
    }

    protected abstract void getG(int Z, int numit, int mult);

    protected abstract void getH(int Z1, int param1, int Z2, int param2);


}