package runcycle;

import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOGeometryOptimization;
import nddoparam.mndo.MNDOSolution;
import nddoparam.param.MNDOParamGradient;
import nddoparam.param.MNDOParamHessian;
import org.jblas.DoubleMatrix;

public abstract class MoleculeRun {
    public String[] output;
    public String hessianStr, newGeomCoords, trainingSet;
    protected MNDOAtom[] atoms, expGeom;
    protected int charge, size, mult;
    protected boolean runHessian;
    protected double[] datum;
    protected MNDOParamGradient g;
    protected MNDOParamHessian h;
    protected MNDOGeometryOptimization opt;
    protected MNDOSolution expSolution;
    protected static double LAMBDA = 1E-7;
    protected String totalHeatDeriv = "";
    protected String totalIonizationDeriv = "";
    protected String totalDipoleDeriv = "";
    protected String totalGeomDeriv = "";
    protected String totalExcelStr = "";
    protected String excelStr = "";
    protected String excelStr2 = "";


    public MoleculeRun(MNDOAtom[] atoms, int charge, MNDOAtom[] expGeom, double[] datum, boolean runHessian, String trainingSet, int mult) {
        this.trainingSet = trainingSet;
        this.atoms = atoms;
        this.expGeom = expGeom;
        this.charge = charge;
        this.runHessian = runHessian;
        this.datum = datum;
        this.size = getSize(trainingSet);
        this.mult = mult;
    }

    private int getSize(String trainingSet) {
        int result = 0;
        int removeH = 0;
        if (trainingSet.contains("H")) {
            result += 5;
            removeH = 1;
        }
        result += (trainingSet.length() - removeH) * 8;
        return result;
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
                    totalHeatDerivSB.append(datum[0]).append(", ").append(g.s.hf).append(", ");
                }


                g.constructErrors(datum[0]);

                if (this.expGeom != null) {
                    g.createExpGeom(this.expGeom, expSolution);
                    g.addGeomError();
                }


                if (datum[1] != 0) {
                    if (numit == 0) {
                        totalDipoleDerivSB.append(datum[1]).append(", ").append(g.s.dipole).append(", ");
                    }

                    g.addDipoleError(datum[1]);
                    totalDipoleDerivSB.append(1 / LAMBDA * (g.eprime.soln.dipole - g.e.soln.dipole)).append(", ");
                }

                if (datum[2] != 0) {
                    if (numit == 0) {
                        totalIonizationDerivSB.append(datum[2]).append(", ").append(-g.s.homo).append(", ");
                    }
                    g.addIEError(datum[1]);
                    totalIonizationDerivSB.append(1 / LAMBDA * (-g.eprime.soln.homo + g.e.soln.homo)).append(", ");
                }

                if (this.expGeom != null) {
                    if (numit == 0) {
                        totalGeomDerivSB.append("0, ").append(g.e.gradient).append(",");
                    }
                    totalGeomDerivSB.append(1 / LAMBDA * (g.eprime.gradient - g.e.gradient)).append(", ");
                }
                
                System.out.println(g.gradient());
                totalHeatDerivSB.append(1 / LAMBDA * (g.eprime.soln.hf - g.e.soln.hf)).append(", ");
                excelSB.append(",").append(g.gradient());
            }
        }
        for (int numit = 0; numit < 8; numit++) {
            getG(6, numit, mult);

            excelStr2 = "," + datum[0] + "," + g.s.hf + ",";


            g.constructErrors(datum[0]);

            if (this.expGeom != null) {
                g.createExpGeom(this.expGeom, expSolution);
                g.addGeomError();
            }

            if (datum[1] != 0) {
                g.addDipoleError(datum[1]);

                excelSB2.append(",").append(datum[1]).append(",").append(g.s.dipole).append(",");

                if (numit == 7 && trainingSet.equals("CH")) {
                    totalDipoleDerivSB.append(1 / LAMBDA * (g.eprime.soln.dipole - g.e.soln.dipole)).append("\n");
                } else {
                    totalDipoleDerivSB.append(1 / LAMBDA * (g.eprime.soln.dipole - g.e.soln.dipole)).append(", ");
                }

            } else {
                excelSB2.append(",,,");
            }

            if (datum[2] != 0) {
                g.addIEError(datum[2]);

                excelSB2.append(",").append(datum[2]).append(",").append(-g.s.homo).append(",");
                if (numit == 7 && trainingSet.equals("CH")) {
                    totalIonizationDerivSB.append(1 / LAMBDA * (-g.eprime.soln.homo + g.e.soln.homo)).append("\n");
                } else {
                    totalIonizationDerivSB.append(1 / LAMBDA * (-g.eprime.soln.homo + g.e.soln.homo)).append(", ");
                }
            } else {
                excelSB2.append(",,,");
            }

            if (this.expGeom != null) {
                excelSB2.append("," + 0 + ",").append(g.e.gradient).append(",");
                if (numit == 7 && trainingSet.equals("CH")) {
                    totalGeomDerivSB.append(1 / LAMBDA * (g.eprime.gradient - g.e.gradient)).append("\n");
                } else {
                    totalGeomDerivSB.append(1 / LAMBDA * (g.eprime.gradient - g.e.gradient)).append(", ");
                }
            }

            if (numit == 7 && trainingSet.equals("CH")) {
                totalHeatDerivSB.append(1 / LAMBDA * (g.eprime.soln.hf - g.e.soln.hf)).append("\n");
            } else {
                totalHeatDerivSB.append(1 / LAMBDA * (g.eprime.soln.hf - g.e.soln.hf)).append(", ");
            }

            excelSB.append(",").append(g.gradient());

        }

        if (trainingSet.strip().contains("N")) {
            for (int numit = 0; numit < 8; numit++) {

                getG(7, numit, mult);


                excelStr2 = "," + datum[0] + "," + g.s.hf + ",";


                g.constructErrors(datum[0]);

                if (this.expGeom != null) {
                    g.createExpGeom(this.expGeom, expSolution);
                    g.addGeomError();
                }

                if (datum[1] != 0) {
                    g.addDipoleError(datum[1]);

                    excelSB2.append(",").append(datum[1]).append(",").append(g.s.dipole).append(",");

                    if (numit == 7) {
                        totalDipoleDerivSB.append(1 / LAMBDA * (g.eprime.soln.dipole - g.e.soln.dipole)).append("\n");
                    } else {
                        totalDipoleDerivSB.append(1 / LAMBDA * (g.eprime.soln.dipole - g.e.soln.dipole)).append(", ");
                    }

                } else {
                    excelSB2.append(",,,");
                }

                if (datum[2] != 0) {
                    g.addIEError(datum[2]);

                    excelSB2.append(",").append(datum[2]).append(",").append(-g.s.homo).append(",");
                    if (numit == 7) {
                        totalIonizationDerivSB.append(1 / LAMBDA * (-g.eprime.soln.homo + g.e.soln.homo)).append("\n");
                    } else {
                        totalIonizationDerivSB.append(1 / LAMBDA * (-g.eprime.soln.homo + g.e.soln.homo)).append(", ");
                    }
                } else {
                    excelSB2.append(",,,");
                }

                if (this.expGeom != null) {
                    excelSB2.append("," + 0 + ",").append(g.e.gradient).append(",");
                    if (numit == 7) {
                        totalGeomDerivSB.append(1 / LAMBDA * (g.eprime.gradient - g.e.gradient)).append("\n");
                    } else {
                        totalGeomDerivSB.append(1 / LAMBDA * (g.eprime.gradient - g.e.gradient)).append(", ");
                    }
                }


                if (numit == 7) {
                    totalHeatDerivSB.append(1 / LAMBDA * (g.eprime.soln.hf - g.e.soln.hf)).append("\n");
                } else {
                    totalHeatDerivSB.append(1 / LAMBDA * (g.eprime.soln.hf - g.e.soln.hf)).append(", ");
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
        System.out.println("Error function: " + g.e.constructErrorFunction());

        totalExcelStr += g.e.constructErrorFunction() + excelStr + excelStr2 + "\n";

        output = new String[]{totalExcelStr, totalHeatDeriv, totalDipoleDeriv, totalIonizationDeriv, totalGeomDeriv};
    }

    protected abstract void getG(int Z, int numit, int mult);

    protected abstract void getH(int Z1, int param1, int Z2, int param2);


}