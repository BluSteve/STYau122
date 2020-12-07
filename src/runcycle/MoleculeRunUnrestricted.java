package runcycle;

import org.jblas.DoubleMatrix;

import mndoparam.mndo.*;
import mndoparam.param.*;

public class MoleculeRunUnrestricted extends AbstractMoleculeRun {

    private MNDOAtom[] atoms;
    private MNDOAtom[] expgeom;
    private MNDOParamGradientUnrestricted g;
    private MNDOParamHessianUnrestricted h;
    private MNDOGeometryOptimizationUnrestricted opt;
    private MNDOSolutionUnrestricted expsoln;


    private static double lambda = 1E-7;

    private static String trainingset = "CHN";

    public MoleculeRunUnrestricted(MNDOAtom[] atms, int charge, int mult, MNDOAtom[] expGeom, double[] datum, boolean runHessian) {
        this.atoms = atms.clone();
        //this.expgeom = expGeom.clone();

        opt = new MNDOGeometryOptimizationUnrestricted(atoms, charge, mult);

        super.newGeomCoords = "UHF\n" + "CHARGE=" + charge + "\nMULT=" + mult + "\n";

        for (int i = 0; i < atoms.length; i++) {
            String s = (i + 1) + "    " + atoms[i].getName() + "    " + String.format("%.9f", atoms[i].getCoordinates()[0] / 1.88973) + "    " + String.format("%.9f", atoms[i].getCoordinates()[1] / 1.88973) + "    " + String.format("%.9f", atoms[i].getCoordinates()[2] / 1.88973) + "\n";
            super.newGeomCoords = super.newGeomCoords + s;

        }

        String totalheatderiv = "";
        String totalionizationderiv = "";
        String totaldipolederiv = "";
        String totalgeomderiv = "";
        String totalexcelstr = "";
        String excelstr = "";
        String excelstr2 = "";

        if (expGeom != null) {
            super.newGeomCoords += "EXPGEOM\n";
            for (int i = 0; i < expGeom.length; i++) {
                String s = (i + 1) + "    " + expGeom[i].getName() + "    " + String.format("%.9f", expGeom[i].getCoordinates()[0] / 1.88973) + "    " + String.format("%.9f", expGeom[i].getCoordinates()[1] / 1.88973) + "    " + String.format("%.9f", expGeom[i].getCoordinates()[2] / 1.88973) + "\n";
                super.newGeomCoords = super.newGeomCoords + s;
            }
            expsoln = new MNDOSolutionUnrestricted(expGeom, charge, mult);
        }

        int size = 13;

        if (trainingset.equals("CHN")) {
            size = 21;
        }


        if (runHessian) {
            DoubleMatrix paramindexes = new DoubleMatrix(size, 2);

            for (int i = 0; i < size; i++) {
                if (i < 5) {
                    paramindexes.put(i, 0, 1);
                } else if (i < 13) {
                    paramindexes.put(i, 0, 6);
                } else if (i < 21) {
                    paramindexes.put(i, 0, 7);
                }
            }

            paramindexes.put(0, 1, 0);
            paramindexes.put(1, 1, 1);
            paramindexes.put(2, 1, 3);
            paramindexes.put(3, 1, 5);
            paramindexes.put(4, 1, 7);

            for (int i = 0; i < 8; i++) {
                paramindexes.put(i + 5, 1, i);
            }
            if (size > 13) {
                for (int i = 0; i < 8; i++) {
                    paramindexes.put(i + 13, 1, i);
                }
            }

            hessianStr = "";

            for (int I = 0; I < size; I++) {
                for (int j = I; j < size; j++) {

                    int Z1 = (int) Math.round(paramindexes.get(I, 0));
                    int Z2 = (int) Math.round(paramindexes.get(j, 0));
                    int param1 = (int) Math.round(paramindexes.get(I, 1));
                    int param2 = (int) Math.round(paramindexes.get(j, 1));

                    h = new MNDOParamHessianUnrestricted(atoms, charge, mult, Z1, param1, Z2, param2, opt.s);


                    h.constructErrors(datum[0]);

                    if (expGeom != null) {
                        h.createExpGeom(expGeom, expsoln);
                        h.addGeomError();
                    }


                    if (datum[1] != 0) {
                        h.addDipoleError(datum[1]);

                    }

                    if (datum[2] != 0) {
                        h.addIEError(datum[2]);

                    }

                    hessianStr += h.hessian() + ", ";

                }
            }

            hessianStr += "\n";


        } else {
            hessianStr = "";
        }


        for (int numit = 0; numit < 8; numit++) {
            if (numit == 2 || numit == 4 || numit == 6) {
            } else {

                g = new MNDOParamGradientUnrestricted(atoms, charge, mult, 1, numit, opt);


                if (numit == 0) {
                    totalheatderiv += datum[0] + ", " + g.s.hf + ", ";
                }


                g.constructErrors(datum[0]);

                if (expGeom != null) {
                    g.createExpGeom(expGeom, expsoln);
                    g.addGeomError();
                }


                if (datum[1] != 0) {
                    if (numit == 0) {
                        totaldipolederiv += datum[1] + ", " + g.s.dipole + ", ";
                    }

                    g.addDipoleError(datum[1]);
                    totaldipolederiv += 1 / lambda * (g.eprime.soln.dipole - g.e.soln.dipole) + ", ";
                }

                if (datum[2] != 0) {
                    if (numit == 0) {
                        totalionizationderiv += datum[2] + ", " + -g.s.homo + ", ";
                    }
                    g.addIEError(datum[1]);
                    totalionizationderiv += 1 / lambda * (-g.eprime.soln.homo + g.e.soln.homo) + ", ";
                }

                if (expGeom != null) {


                    if (numit == 0) {
                        totalgeomderiv += "0, " + g.e.gradient + ",";
                    }
                    totalgeomderiv += 1 / lambda * (g.eprime.gradient - g.e.gradient) + ", ";
                }


                System.out.println(g.gradient());


                totalheatderiv += 1 / lambda * (g.eprime.soln.hf - g.e.soln.hf) + ", ";


                excelstr += "," + g.gradient();


            }

        }

        for (int numit = 0; numit < 8; numit++) {

            g = new MNDOParamGradientUnrestricted(atoms, charge, mult, 6, numit, opt);


            excelstr2 = "," + datum[0] + "," + g.s.hf + ",";


            g.constructErrors(datum[0]);

            if (expGeom != null) {
                g.createExpGeom(expGeom, expsoln);
                g.addGeomError();
            }

            if (datum[1] != 0) {
                g.addDipoleError(datum[1]);

                excelstr2 += "," + datum[1] + "," + g.s.dipole + ",";

                if (numit == 7 && trainingset.equals("CH")) {
                    totaldipolederiv += 1 / lambda * (g.eprime.soln.dipole - g.e.soln.dipole) + "\n";
                } else {
                    totaldipolederiv += 1 / lambda * (g.eprime.soln.dipole - g.e.soln.dipole) + ", ";
                }

            } else {
                excelstr2 += ",,,";
            }

            if (datum[2] != 0) {
                g.addIEError(datum[2]);

                excelstr2 += "," + datum[2] + "," + -g.s.homo + ",";
                if (numit == 7 && trainingset.equals("CH")) {
                    totalionizationderiv += 1 / lambda * (-g.eprime.soln.homo + g.e.soln.homo) + "\n";
                } else {
                    totalionizationderiv += 1 / lambda * (-g.eprime.soln.homo + g.e.soln.homo) + ", ";
                }
            } else {
                excelstr2 += ",,,";
            }

            if (expGeom != null) {
                excelstr2 += "," + 0 + "," + g.e.gradient + ",";
                if (numit == 7 && trainingset.equals("CH")) {
                    totalgeomderiv += 1 / lambda * (g.eprime.gradient - g.e.gradient) + "\n";
                } else {
                    totalgeomderiv += 1 / lambda * (g.eprime.gradient - g.e.gradient) + ", ";
                }
            }


            if (numit == 7 && trainingset.equals("CH")) {
                totalheatderiv += 1 / lambda * (g.eprime.soln.hf - g.e.soln.hf) + "\n";
            } else {
                totalheatderiv += 1 / lambda * (g.eprime.soln.hf - g.e.soln.hf) + ", ";
            }

            excelstr += "," + g.gradient();

        }

        if (trainingset.strip().contains("N")) {
            for (int numit = 0; numit < 8; numit++) {

                g = new MNDOParamGradientUnrestricted(atoms, charge, mult, 7, numit, opt);


                excelstr2 = "," + datum[0] + "," + g.s.hf + ",";


                g.constructErrors(datum[0]);

                if (expGeom != null) {
                    g.createExpGeom(expGeom, expsoln);
                    g.addGeomError();
                }

                if (datum[1] != 0) {
                    g.addDipoleError(datum[1]);

                    excelstr2 += "," + datum[1] + "," + g.s.dipole + ",";

                    if (numit == 7) {
                        totaldipolederiv += 1 / lambda * (g.eprime.soln.dipole - g.e.soln.dipole) + "\n";
                    } else {
                        totaldipolederiv += 1 / lambda * (g.eprime.soln.dipole - g.e.soln.dipole) + ", ";
                    }

                } else {
                    excelstr2 += ",,,";
                }

                if (datum[2] != 0) {
                    g.addIEError(datum[2]);

                    excelstr2 += "," + datum[2] + "," + -g.s.homo + ",";
                    if (numit == 7) {
                        totalionizationderiv += 1 / lambda * (-g.eprime.soln.homo + g.e.soln.homo) + "\n";
                    } else {
                        totalionizationderiv += 1 / lambda * (-g.eprime.soln.homo + g.e.soln.homo) + ", ";
                    }
                } else {
                    excelstr2 += ",,,";
                }

                if (expGeom != null) {
                    excelstr2 += "," + 0 + "," + g.e.gradient + ",";
                    if (numit == 7) {
                        totalgeomderiv += 1 / lambda * (g.eprime.gradient - g.e.gradient) + "\n";
                    } else {
                        totalgeomderiv += 1 / lambda * (g.eprime.gradient - g.e.gradient) + ", ";
                    }
                }


                if (numit == 7) {
                    totalheatderiv += 1 / lambda * (g.eprime.soln.hf - g.e.soln.hf) + "\n";
                } else {
                    totalheatderiv += 1 / lambda * (g.eprime.soln.hf - g.e.soln.hf) + ", ";
                }

                excelstr += "," + g.gradient();

            }
        }


        System.out.println("Error function: " + g.e.constructErrorFunction());

        totalexcelstr += g.e.constructErrorFunction() + excelstr + excelstr2 + "\n";

        output = new String[]{totalexcelstr, totalheatderiv, totaldipolederiv, totalionizationderiv, totalgeomderiv};

    }

}
