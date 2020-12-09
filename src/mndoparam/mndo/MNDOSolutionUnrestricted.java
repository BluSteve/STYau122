package mndoparam.mndo;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;

import java.util.ArrayList;
import java.util.Arrays;


public class MNDOSolutionUnrestricted extends AbstractMNDOSolution {
    private DoubleMatrix J, Ka, Kb, Fa, Fb, Ca, Cb;
    DoubleMatrix Ea;
    DoubleMatrix Eb;

    private DoubleMatrix alphaDensity;
    private DoubleMatrix betaDensity;
    private int Nalpha, Nbeta;

    private double[] integralArrayCoulomb, integralArrayExchange;

    public int multiplicity;

    public MNDOSolutionUnrestricted(MNDOAtom[] atoms, int charge, int mult) {
        super(atoms, charge);
        this.multiplicity = mult;
        if (nElectrons % 2 == multiplicity % 2 || multiplicity < 1) {
            System.err.println("You're high. (Please check multiplicity and charge): " + nElectrons + ", " + multiplicity);
            System.exit(0);
        }

        nElectrons -= (multiplicity - 1);

        Nalpha = nElectrons / 2 + (multiplicity - 1);

        Nbeta = nElectrons / 2;

        System.out.println("1-electron matrix elements evaluated - moving on to two-electron matrix");


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
                }
            }
        }
        integralArrayCoulomb = new double[size];
        int integralCount = 0;
        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                if (j == k) {
                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            integralArrayCoulomb[integralCount] = MNDO6G.OneCenterERI(orbitals[j], orbitals[j], orbitals[l], orbitals[l]);
                            integralCount++;
                        }
                    }
                    for (int l : missingIndex[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : missingIndex[atomNumber[j]]) {
                                if (m > -1) {
                                    if (atomNumber[l] == atomNumber[m]) {
                                        integralArrayCoulomb[integralCount] = MNDO6G.getG(orbitals[j], orbitals[j], orbitals[l], orbitals[m]);
                                        integralCount++;
                                    }
                                }
                            }
                        }
                    }
                } else if (atomNumber[j] == atomNumber[k]) {
                    integralArrayCoulomb[integralCount] = 2 * MNDO6G.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]);
                    integralCount++;
                    for (int l : missingIndex[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : missingIndex[atomNumber[j]]) {
                                if (m > -1) {
                                    if (atomNumber[l] == atomNumber[m]) {
                                        integralArrayCoulomb[integralCount] = MNDO6G.getG(orbitals[j], orbitals[k], orbitals[l], orbitals[m]);
                                        integralCount++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        System.out.println("Coulomb (J) matrix ERIs evaluated - moving on to Exchange (K) matrix ERIs...");

        size = 0;
        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                //System.err.println ("(" + j + ", " + k + ")");
                if (j == k) {

                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            size++;
                        }
                    }
                } else if (atomNumber[j] == atomNumber[k]) {
                    size++;
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

        integralArrayExchange = new double[size];
        integralCount = 0;
        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                //System.err.println ("(" + j + ", " + k + ")");
                if (j == k) {

                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            integralArrayExchange[integralCount] = -1 * MNDO6G.OneCenterERI(orbitals[j], orbitals[l], orbitals[j], orbitals[l]);
                            integralCount++;
                        }
                    }
                } else if (atomNumber[j] == atomNumber[k]) {
                    //System.err.println ("1.5[" + j + k + "|" + j + k + "] - 0.5[" + j + j + "|" + k + k + "]");
                    integralArrayExchange[integralCount] = -1 * MNDO6G.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]) - 1 * MNDO6G.OneCenterERI(orbitals[j], orbitals[j], orbitals[k], orbitals[k]);
                    integralCount++;
                } else {
                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : index[atomNumber[k]]) {
                                if (m > -1) {
                                    integralArrayExchange[integralCount] = -1 * MNDO6G.getG(orbitals[j], orbitals[l], orbitals[k], orbitals[m]);
                                    integralCount++;
                                }
                            }
                        }
                    }
                }
            }
        }

        DoubleMatrix[] matrices = Eigen.symmetricEigenvectors(H);

        System.out.println("Exchange (K) matrix ERIs evaluated, beginning SCF iterations...");

        Ea = matrices[1].diag();

        Eb = matrices[1].diag();

        Ca = matrices[0].transpose();

        Cb = matrices[0].transpose();

        J = new DoubleMatrix(Ca.rows, Ca.columns);

        Ka = new DoubleMatrix(Ca.rows, Ca.columns);

        Kb = new DoubleMatrix(Ca.rows, Ca.columns);

        alphaDensity = DensityMatrix(Ca, Nalpha);

        betaDensity = DensityMatrix(Cb, Nbeta);

        DoubleMatrix oldalphadensity = DoubleMatrix.zeros(Ca.rows, Ca.columns);

        DoubleMatrix oldbetadensity = DoubleMatrix.zeros(Ca.rows, Ca.columns);

        int Jcount, Kcount;

        int numIt = 0;
        while (!(isSimilar(alphaDensity, oldalphadensity, 1E-10) && isSimilar(betaDensity, oldbetadensity, 1E-10))) {

            numIt++;
            oldalphadensity = alphaDensity.dup();
            oldbetadensity = betaDensity.dup();

            Jcount = 0;
            Kcount = 0;

            //construct J matrix

            for (int j = 0; j < orbitals.length; j++) {
                for (int k = j; k < orbitals.length; k++) {
                    double val = 0;
                    if (j == k) {

                        for (int l : index[atomNumber[j]]) {
                            if (l > -1) {
                                val += (alphaDensity.get(l, l) + betaDensity.get(l, l)) * integralArrayCoulomb[Jcount];
                                Jcount++;
                            }
                        }

                        for (int l : missingIndex[atomNumber[j]]) {
                            if (l > -1) {
                                for (int m : missingIndex[atomNumber[j]]) {
                                    if (m > -1) {
                                        if (atomNumber[l] == atomNumber[m]) {
                                            val += (alphaDensity.get(l, m) + betaDensity.get(l, m)) * integralArrayCoulomb[Jcount];
                                            Jcount++;
                                        }
                                    }

                                }
                            }
                        }
                    } else if (atomNumber[j] == atomNumber[k]) {
                        val += (alphaDensity.get(j, k) + betaDensity.get(j, k)) * integralArrayCoulomb[Jcount];
                        Jcount++;

                        for (int l : missingIndex[atomNumber[j]]) {
                            if (l > -1) {
                                for (int m : missingIndex[atomNumber[j]]) {
                                    if (m > -1) {
                                        if (atomNumber[l] == atomNumber[m]) {
                                            val += (alphaDensity.get(l, m) + betaDensity.get(l, m)) * integralArrayCoulomb[Jcount];
                                            Jcount++;
                                        }
                                    }

                                }
                            }
                        }
                    }


                    J.put(j, k, val);
                    J.put(k, j, val);
                }
            }

            for (int j = 0; j < orbitals.length; j++) {
                for (int k = j; k < orbitals.length; k++) {
                    double vala = 0;
                    double valb = 0;
                    if (j == k) {

                        for (int l : index[atomNumber[j]]) {
                            if (l > -1) {
                                vala += alphaDensity.get(l, l) * integralArrayExchange[Kcount];
                                valb += betaDensity.get(l, l) * integralArrayExchange[Kcount];
                                Kcount++;
                            }
                        }

                    } else if (atomNumber[j] == atomNumber[k]) {
                        vala += alphaDensity.get(j, k) * integralArrayExchange[Kcount];
                        valb += betaDensity.get(j, k) * integralArrayExchange[Kcount];
                        Kcount++;

                    } else {
                        for (int l : index[atomNumber[j]]) {
                            if (l > -1) {
                                for (int m : index[atomNumber[k]]) {
                                    if (m > -1) {
                                        vala += alphaDensity.get(l, m) * integralArrayExchange[Kcount];
                                        valb += betaDensity.get(l, m) * integralArrayExchange[Kcount];
                                        Kcount++;
                                    }
                                }
                            }
                        }
                    }

                    Ka.put(j, k, vala);
                    Ka.put(k, j, vala);
                    Kb.put(j, k, valb);
                    Kb.put(k, j, valb);
                }
            }

            Fa = H.add(J).add(Ka);
            Fb = H.add(J).add(Kb);

            DoubleMatrix[] matrices1 = Eigen.symmetricEigenvectors(Fa);

            DoubleMatrix[] matrices2 = Eigen.symmetricEigenvectors(Fb);

            Ea = matrices1[1].diag();

            Eb = matrices2[1].diag();

            Ca = matrices1[0].transpose();

            Cb = matrices2[0].transpose();

            alphaDensity = DensityMatrix(Ca, Nalpha).mmul(1 - damp).add(oldalphadensity.mmul(damp));

            betaDensity = DensityMatrix(Cb, Nbeta).mmul(1 - damp).add(oldbetadensity.mmul(damp));

            if (numIt >= 100000) {
                System.err.println("SCF Has Not Converged");

                System.err.println("Damping Coefficient will be Increased, and the run restarted...");

                damp += 0.02;

                matrices = Eigen.symmetricEigenvectors(H);

                System.out.println("Exchange (K) matrix ERIs evaluated, beginning SCF iterations...");

                Ea = matrices[1].diag();

                Eb = matrices[1].diag();

                Ca = matrices[0].transpose();

                Cb = matrices[0].transpose();

                J = new DoubleMatrix(Ca.rows, Ca.columns);

                Ka = new DoubleMatrix(Ca.rows, Ca.columns);

                Kb = new DoubleMatrix(Ca.rows, Ca.columns);

                alphaDensity = DensityMatrix(Ca, Nalpha);

                betaDensity = DensityMatrix(Cb, Nbeta);

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

        //System.out.println ("Alpha eigenvalues: " + Ea);
        //System.out.println ("Beta eigenvalues: " + Eb);

        double e = 0;

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = 0; k < orbitals.length; k++) {
                e += 0.5 * alphaDensity.get(j, k) * (H.get(j, k) + Fa.get(j, k));
                e += 0.5 * betaDensity.get(j, k) * (H.get(j, k) + Fb.get(j, k));
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

        this.homo = Ea.get(Nalpha - 1, 0);
        this.lumo = 0.001 * Math.round(Eb.get(Nbeta, 0) * 1000);

        //System.out.println ("HOMO energy: " + homo + " eV");

        //System.out.println ("LUMO energy: " + lumo + " eV");

        //System.out.println ("Heat of Formation: " + 0.01 * Math.round(heat * 100) + " eV = " + 0.01 * Math.round(heat * 100 / 4.3363E-2) + "kcal/mol");

        double[] populations = new double[atoms.length];

        for (int j = 0; j < atoms.length; j++) {
            double sum = 0;
            for (int k : index[j]) {
                if (k > -1) {
                    sum += alphaDensity.get(k, k) + betaDensity.get(k, k);
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

        //chargedip[0] = 0.001 * Math.round(chargedip[0] * 1000);
        //chargedip[1] = 0.001 * Math.round(chargedip[1] * 1000);
        //chargedip[2] = 0.001 * Math.round(chargedip[2] * 1000);

        //System.out.println ("point charge dipole contribution: " + Arrays.toString(chargedip));

        hybridip = new double[]{0, 0, 0};

        for (int j = 0; j < atoms.length; j++) {

            if (index[j][1] != -1) {//exclude hydrogen
                hybridip[0] = hybridip[0] - 2.5416 * 2 * atoms[j].D1 * (alphaDensity.get(index[j][0], index[j][1]) + betaDensity.get(index[j][0], index[j][1]));
                hybridip[1] = hybridip[1] - 2.5416 * 2 * atoms[j].D1 * (alphaDensity.get(index[j][0], index[j][2]) + betaDensity.get(index[j][0], index[j][2]));
                hybridip[2] = hybridip[2] - 2.5416 * 2 * atoms[j].D1 * (alphaDensity.get(index[j][0], index[j][3]) + betaDensity.get(index[j][0], index[j][3]));
            }
        }

        //hybridip[0] = 0.001 * Math.round(hybridip[0] * 1000);
        //hybridip[1] = 0.001 * Math.round(hybridip[1] * 1000);
        //hybridip[2] = 0.001 * Math.round(hybridip[2] * 1000);

        //System.out.println ("hybrid dipole contribution: " + Arrays.toString(hybridip));

        dipoletot = new double[]{chargedip[0] + hybridip[0], chargedip[1] + hybridip[1], chargedip[2] + hybridip[2]};

        //System.out.println ("sum: " + Arrays.toString(dipoletot));


        dipole = Math.sqrt(dipoletot[0] * dipoletot[0] + dipoletot[1] * dipoletot[1] + dipoletot[2] * dipoletot[2]);

        //System.out.println ("dipole moment: " + dipole + " debyes");
    }

    private DoubleMatrix DensityMatrix(DoubleMatrix c, int NElectrons) {//density matrix construction by definition.


        DoubleMatrix densitymatrix = new DoubleMatrix(orbitals.length, orbitals.length);

        for (int i = 0; i < orbitals.length; i++) {
            for (int j = 0; j < orbitals.length; j++) {
                double sum = 0;
                int count = NElectrons;
                int counter = -1;


                while (count > 0) {

                    counter++;

                    sum += c.get(counter, i) * c.get(counter, j);
                    count -= 1;

                }


                densitymatrix.put(i, j, sum);
            }
        }


        return densitymatrix;
    }

    public DoubleMatrix alphadensity() {
        return this.alphaDensity;
    }

    public DoubleMatrix betadensity() {
        return this.betaDensity;
    }
}
