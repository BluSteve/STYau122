package mndoparam.mndo;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;

import java.util.Arrays;


public class MNDOSolution extends AbstractMNDOSolution {

    private DoubleMatrix densityMatrix;

    public DoubleMatrix C, F, G, E;//H - core matrix, G = 2-electron matrix, F = fock matrix, C = coeffecient matrix (transposed for easier reading), E = eigenvalues
    public MNDOSolution(MNDOAtom[] atoms, int charge) {
        super(atoms, charge);
        int size = 0;
        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                if (j == k) {

                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            size++;
                        }
                    }

                    for (int l : missingIndex[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : missingIndex[atomNumber[j]]) {
                                if (m > -1) {
                                    if (atomNumber[l] == atomNumber[m]) {
                                        size++;
                                    }
                                }

                            }
                        }
                    }
                } else if (atomNumber[j] == atomNumber[k]) {
                    size++;

                    for (int l : missingIndex[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : missingIndex[atomNumber[j]]) {
                                if (m > -1) {
                                    if (atomNumber[l] == atomNumber[m]) {
                                        size++;
                                    }
                                }

                            }
                        }
                    }
                } else {
                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : index[atomNumber[k]]) {
                                if (m > -1) {
                                    size++;
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.println("1-electron matrix elements evaluated - moving on to two-electron matrix");
        double[] integralArray = new double[size];
        //The idea of the integralarray is to simply store all the integrals in order they are called. It's basically my way of avoiding having to perform a Yoshemine sort.

        int integralcount = 0;
        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                if (j == k) {

                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            integralArray[integralcount] = (MNDO6G.OneCenterERI(orbitals[j], orbitals[j], orbitals[l], orbitals[l]) - 0.5 * MNDO6G.OneCenterERI(orbitals[j], orbitals[l], orbitals[j], orbitals[l]));
                            integralcount++;
                        }
                    }

                    for (int l : missingIndex[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : missingIndex[atomNumber[j]]) {
                                if (m > -1) {
                                    if (atomNumber[l] == atomNumber[m]) {
                                        integralArray[integralcount] = (MNDO6G.getG(orbitals[j], orbitals[j], orbitals[l], orbitals[m]));
                                        integralcount++;
                                    }
                                }

                            }
                        }
                    }
                } else if (atomNumber[j] == atomNumber[k]) {
                    integralArray[integralcount] = (1.5 * MNDO6G.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]) - 0.5 * MNDO6G.OneCenterERI(orbitals[j], orbitals[j], orbitals[k], orbitals[k]));
                    integralcount++;
                    for (int l : missingIndex[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : missingIndex[atomNumber[j]]) {
                                if (m > -1) {
                                    if (atomNumber[l] == atomNumber[m]) {
                                        integralArray[integralcount] = (MNDO6G.getG(orbitals[j], orbitals[k], orbitals[l], orbitals[m]));
                                        integralcount++;
                                    }
                                }

                            }
                        }
                    }
                } else {
                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : index[atomNumber[k]]) {
                                if (m > -1) {
                                    integralArray[integralcount] = (-0.5 * MNDO6G.getG(orbitals[j], orbitals[l], orbitals[k], orbitals[m]));
                                    integralcount++;
                                }
                            }
                        }
                    }
                }
            }
        }
        System.out.println("2-electron integrals evaluated");

        DoubleMatrix[] matrices = Eigen.symmetricEigenvectors(H);

        System.out.println("initial diagonalization completed, beginning SCF iterations...");

        E = matrices[1].diag();

        C = matrices[0].transpose();

        G = DoubleMatrix.zeros(C.rows, C.columns);

        densityMatrix = DensityMatrix(C);

        DoubleMatrix olddensity = DoubleMatrix.zeros(C.rows, C.columns);

        F = H.dup();


        int numIt = 0;

        while (!isSimilar(densityMatrix, olddensity, 1E-10)) {//density matrix convergence criteria; since each iteration takes place within a fraction of a second I figured why not

            numIt++;
            olddensity = densityMatrix.dup();

            integralcount = 0;

            //this entire block of code fills up the G matrix, and it calls the integralarray to save time.

            for (int j = 0; j < orbitals.length; j++) {
                for (int k = j; k < orbitals.length; k++) {
                    double val = 0;
                    if (j == k) {

                        for (int l : index[atomNumber[j]]) {
                            if (l > -1) {
                                val += densityMatrix.get(l, l) * integralArray[integralcount];
                                integralcount++;
                            }
                        }

                        for (int l : missingIndex[atomNumber[j]]) {
                            if (l > -1) {
                                for (int m : missingIndex[atomNumber[j]]) {
                                    if (m > -1) {
                                        if (atomNumber[l] == atomNumber[m]) {
                                            val += densityMatrix.get(l, m) * integralArray[integralcount];
                                            integralcount++;
                                        }
                                    }

                                }
                            }
                        }
                    } else if (atomNumber[j] == atomNumber[k]) {
                        val += densityMatrix.get(j, k) * integralArray[integralcount];
                        integralcount++;

                        for (int l : missingIndex[atomNumber[j]]) {
                            if (l > -1) {
                                for (int m : missingIndex[atomNumber[j]]) {
                                    if (m > -1) {
                                        if (atomNumber[l] == atomNumber[m]) {
                                            val += densityMatrix.get(l, m) * integralArray[integralcount];
                                            integralcount++;
                                        }
                                    }

                                }
                            }
                        }
                    } else {
                        for (int l : index[atomNumber[j]]) {
                            if (l > -1) {
                                for (int m : index[atomNumber[k]]) {
                                    if (m > -1) {
                                        val += densityMatrix.get(l, m) * integralArray[integralcount];
                                        integralcount++;
                                    }
                                }
                            }
                        }
                    }

                    G.put(j, k, val);
                    G.put(k, j, val);
                }
            }

            F = H.dup().add(G);


            matrices = Eigen.symmetricEigenvectors(F);

            E = matrices[1].diag();

            //System.err.println (E);

            C = matrices[0].transpose();

            densityMatrix = DensityMatrix(C).mmul(1 - damp).add(olddensity.mmul(damp));

            //System.out.println (densitymatrix);

            if (numIt >= 10000) {
                System.err.println("SCF Has Not Converged");

                System.err.println("Damping Coefficient will be Increased, and the run restarted...");

                damp += 0.02;

                matrices = Eigen.symmetricEigenvectors(H);


                E = matrices[1].diag();

                C = matrices[0].transpose();

                G = DoubleMatrix.zeros(C.rows, C.columns);

                densityMatrix = DensityMatrix(C);

                numIt = 0;

                if (damp >= 1) {
                    System.err.println("Damping Coefficient Cannot Be Increased Further. Exiting program...");

                    for (MNDOAtom a : atoms) {
                        System.out.println(a.getAtomProperties().getZ() + "; " + Arrays.toString(a.getCoordinates()));
                    }
                    System.exit(0);

                }
            }

        }

        System.out.println("SCF completed");

        //System.out.println ("Eigenvalues: " + E + "\n\n");

        //System.out.println ("Ionization energy: " + 0.001 * Math.round(E.get(nElectrons / 2 - 1, 0) * 1000) + " eV");

        //System.out.println ("LUMO energy: " + 0.001 * Math.round(E.get(nElectrons / 2, 0) * 1000) + " eV");

        double e = 0;

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = 0; k < orbitals.length; k++) {
                e += 0.5 * densityMatrix.get(j, k) * (H.get(j, k) + F.get(j, k));
            }
        }
        //System.out.println ("Electronic energy: " + 0.01 * Math.round(e * 100) + " eV");
        double heat = 0;

        for (int j = 0; j < atoms.length; j++) {
            heat += atoms[j].getHeat() - atoms[j].getEisol();
            for (int k = j + 1; k < atoms.length; k++) {
                e += MNDOAtom.crf(atoms[j], atoms[k]);
            }
        }
        //System.out.println ("Core repulsion energy: " + 0.01 * Math.round(d * 100) + " eV");
        //System.out.println ("Energy: " + 0.01 * Math.round(e * 100) + " eV");

        energy = e;

        heat += e;

        this.hf = heat / 4.3363E-2;

        if (nElectrons > 0) {
            this.homo = E.get(nElectrons / 2 - 1, 0);
        } else {
            this.homo = 0;
        }

        this.lumo = E.get(nElectrons / 2, 0);

        //System.out.println ("Heat of Formation: " + 0.01 * Math.round(heat * 100) + " eV = " + 0.01 * Math.round(heat * 100 / 4.3363E-2) + "kcal/mol");

        double[] populations = new double[atoms.length];

        for (int j = 0; j < atoms.length; j++) {
            double sum = 0;
            for (int k : index[j]) {
                if (k > -1) {
                    sum += densityMatrix.get(k, k);
                }
            }

            populations[j] = atoms[j].getAtomProperties().getQ() - sum;
        }


        double[] com = new double[]{0, 0, 0};

        double mass = 0;

        for (MNDOAtom atom : atoms) {
            com[0] = com[0] + atom.getMass() * atom.getCoordinates()[0];
            com[1] = com[1] + atom.getMass() * atom.getCoordinates()[1];
            com[2] = com[2] + atom.getMass() * atom.getCoordinates()[2];
            mass += atom.getMass();
        }

        com[0] = com[0] / mass;
        com[1] = com[1] / mass;
        com[2] = com[2] / mass;


        chargedip = new double[]{0, 0, 0};

        for (int j = 0; j < atoms.length; j++) {
            chargedip[0] += 2.5416 * populations[j] * (atoms[j].getCoordinates()[0] - com[0]);
            chargedip[1] += 2.5416 * populations[j] * (atoms[j].getCoordinates()[1] - com[1]);
            chargedip[2] += 2.5416 * populations[j] * (atoms[j].getCoordinates()[2] - com[2]);
        }

        //chargedip[0] = 0.01 * Math.round(chargedip[0] * 100);
        //chargedip[1] = 0.01 * Math.round(chargedip[1] * 100);
        //chargedip[2] = 0.01 * Math.round(chargedip[2] * 100);

        //System.out.println ("point charge dipole contribution: " + Arrays.toString(chargedip));

        hybridip = new double[]{0, 0, 0};

        for (int j = 0; j < atoms.length; j++) {

            if (index[j][1] != -1) {//exclude hydrogen
                hybridip[0] = hybridip[0] - 2.5416 * 2 * atoms[j].D1 * densityMatrix.get(index[j][0], index[j][1]);
                hybridip[1] = hybridip[1] - 2.5416 * 2 * atoms[j].D1 * densityMatrix.get(index[j][0], index[j][2]);
                hybridip[2] = hybridip[2] - 2.5416 * 2 * atoms[j].D1 * densityMatrix.get(index[j][0], index[j][3]);
            }
        }

        //hybridip[0] = 0.01 * Math.round(hybridip[0] * 100);
        //hybridip[1] = 0.01 * Math.round(hybridip[1] * 100);
        //hybridip[2] = 0.01 * Math.round(hybridip[2] * 100);

        //System.out.println ("hybrid dipole contribution: " + Arrays.toString(hybridip));

        dipoletot = new double[]{chargedip[0] + hybridip[0], chargedip[1] + hybridip[1], chargedip[2] + hybridip[2]};

        //System.out.println ("sum: " + Arrays.toString(dipoletot));


        dipole = Math.sqrt(dipoletot[0] * dipoletot[0] + dipoletot[1] * dipoletot[1] + dipoletot[2] * dipoletot[2]);

        //System.out.println ("dipole moment: " + dipole + " debyes");

    }

    public DoubleMatrix getE() {
        return E;
    }


    private DoubleMatrix DensityMatrix(DoubleMatrix c) {//density matrix construction by definition.


        DoubleMatrix densitymatrix = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int i = 0; i < orbitals.length; i++) {
            for (int j = 0; j < orbitals.length; j++) {
                double sum = 0;
                int count = nElectrons;
                int counter = -1;


                while (count > 0) {

                    counter++;

                    sum += 2 * c.get(counter, i) * c.get(counter, j);
                    count -= 2;

                }


                densitymatrix.put(i, j, sum);
            }
        }
        return densitymatrix;
    }



    public DoubleMatrix densitymatrix() {
        return this.densityMatrix;
    }
}
