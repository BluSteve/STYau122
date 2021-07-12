package nddoparam;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.jblas.Solve;
import org.apache.commons.lang3.time.StopWatch;

import java.util.Arrays;


public class NDDOSolutionRestrictedDIIS extends NDDOSolution {

    private DoubleMatrix densityMatrix;

    public double[] integralArray;

    public DoubleMatrix C, F, G, E;//H - core matrix, G = 2-electron matrix, F = fock matrix, C = coeffecient matrix (transposed for easier reading), E = eigenvalues

    public NDDOSolutionRestrictedDIIS(NDDOAtom[] atoms, int charge) {
        super(atoms, charge);

        StopWatch sw = new StopWatch();

        sw.start();

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
        integralArray = new double[size];
        //The idea of the integralarray is to simply store all the integrals in order they are called. It's basically my way of avoiding having to perform a Yoshemine sort.
        // TODO re-implement HashMap
        int integralcount = 0;
        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                if (j == k) { // case 1
                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            integralArray[integralcount] = (NDDO6G.OneCenterERI(orbitals[j], orbitals[j], orbitals[l], orbitals[l]) - 0.5 * NDDO6G.OneCenterERI(orbitals[j], orbitals[l], orbitals[j], orbitals[l]));
                            integralcount++;
                        }
                    }

                    for (int l : missingIndex[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : missingIndex[atomNumber[j]]) {
                                if (m > -1) {
                                    if (atomNumber[l] == atomNumber[m]) {
                                        integralArray[integralcount] = (NDDO6G.getG(orbitals[j], orbitals[j], orbitals[l], orbitals[m]));
                                        integralcount++;
                                    }
                                }

                            }
                        }
                    }
                } else if (atomNumber[j] == atomNumber[k]) { // case 2
                    integralArray[integralcount] = (1.5 * NDDO6G.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]) - 0.5 * NDDO6G.OneCenterERI(orbitals[j], orbitals[j], orbitals[k], orbitals[k]));
                    integralcount++;
                    for (int l : missingIndex[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : missingIndex[atomNumber[j]]) {
                                if (m > -1) {
                                    if (atomNumber[l] == atomNumber[m]) {
                                        integralArray[integralcount] = (NDDO6G.getG(orbitals[j], orbitals[k], orbitals[l], orbitals[m]));
                                        integralcount++;
                                    }
                                }
                            }
                        }
                    }
                } else { // case 3
                    for (int l : index[atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : index[atomNumber[k]]) {
                                if (m > -1) {
                                    integralArray[integralcount] = (-0.5 * NDDO6G.getG(orbitals[j], orbitals[l], orbitals[k], orbitals[m]));
                                    integralcount++;
                                }
                            }
                        }
                    }
                }
            }
        }

        DoubleMatrix[] matrices = Eigen.symmetricEigenvectors(H);

        System.out.println(moleculeName + " Initial diagonalization completed, beginning SCF iterations...");

        E = matrices[1].diag();

        C = matrices[0].transpose();

        G = DoubleMatrix.zeros(C.rows, C.columns);

        densityMatrix = calculateDensityMatrix(C);

        DoubleMatrix olddensity = DoubleMatrix.zeros(C.rows, C.columns);

        F = H.dup();

        DoubleMatrix[] Farray = new DoubleMatrix[8];
        DoubleMatrix[] Darray = new DoubleMatrix[8];
        DoubleMatrix[] commutatorarray = new DoubleMatrix[8];

        DoubleMatrix B = DoubleMatrix.zeros (8, 8);


        int numIt = 0;

        double DIISError = 10;
        while (DIISError > 1E-10) {//density matrix convergence criteria; since each iteration takes place within a fraction of a second I figured why not



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

            if (numIt < Farray.length) {

                Farray[numIt] = F.dup();

                Darray[numIt] = densityMatrix.dup();
                commutatorarray[numIt] = commutator(F.dup(), densityMatrix.dup());
                DIISError = commutatorarray[numIt].norm2();

                for (int i = 0; i <= numIt; i++) {

                    double product = (commutatorarray[numIt].mmul(commutatorarray[i].transpose())).diag().sum();
                    B.put(i, numIt, product);
                    B.put(numIt, i, product);
                }
            }
            else {

                for (int i = 0; i < Farray.length - 1; i++) {

                    Farray[i] = Farray[i+1].dup();
                    Darray[i] = Darray[i+1].dup();
                    commutatorarray[i] = commutatorarray[i+1].dup();
                }

                Farray[Farray.length - 1] = F.dup();
                Darray[Darray.length - 1] = densityMatrix.dup();
                commutatorarray[Darray.length - 1] = commutator(F.dup(), densityMatrix.dup());
                DIISError = commutatorarray[Darray.length - 1].norm2();

                DoubleMatrix newB = DoubleMatrix.zeros (8, 8);

                for (int i = 0; i < Farray.length - 1; i++) {
                    for (int j = i; j < Farray.length - 1; j++) {
                        newB.put(i, j, B.get(i + 1, j + 1));
                        newB.put(j, i , B.get(i + 1, j + 1));
                    }
                }

                for (int i = 0; i < Farray.length; i++) {

                    double product = commutatorarray[Farray.length - 1].transpose().mmul(commutatorarray[i]).diag().sum();
                    newB.put(i, Farray.length - 1, product);
                    newB.put(Farray.length - 1, i, product);
                }

                B = newB.dup();
            }


            if (commutatorarray[Math.min(Farray.length - 1, numIt)].max() < 0.1) {

                DoubleMatrix mat = DoubleMatrix.zeros (Math.min(Farray.length + 1, numIt + 2), Math.min(Farray.length + 1, numIt + 2));

                for (int i = 0; i < Math.min(Farray.length, numIt + 1); i++) {
                    for (int j = i; j < Math.min(Farray.length, numIt + 1); j++) {
                        mat.put(i, j, B.get(i, j));
                        mat.put(j, i, B.get(i, j));

                    }
                }



                mat.putColumn(mat.columns - 1, DoubleMatrix.ones(mat.rows, 1));

                mat.putRow(mat.rows - 1, DoubleMatrix.ones(mat.columns, 1));

                mat.put(mat.rows - 1, mat.columns - 1, 0);

                DoubleMatrix rhs = DoubleMatrix.zeros (mat.rows, 1);

                rhs.put(mat.rows - 1, 0, 1);

                try {
                    DoubleMatrix DIIS = Solve.solve(mat, rhs);

                    DoubleMatrix F = DoubleMatrix.zeros (densityMatrix.rows, densityMatrix.columns);



                    for (int i = 0; i < DIIS.length - 1; i++) {
                        F = F.add (Farray[i].mmul(DIIS.get(i)));
                    }


                    this.F = F.dup();


                    matrices = Eigen.symmetricEigenvectors(F);

                    E = matrices[1].diag();

                    C = matrices[0].transpose();

                    densityMatrix = calculateDensityMatrix(C);
                }
                catch (Exception e) {
                    matrices = Eigen.symmetricEigenvectors(F);

                    E = matrices[1].diag();

                    C = matrices[0].transpose();

                    densityMatrix = calculateDensityMatrix(C).mmul(1 - damp).add(olddensity.mmul(damp));
                }


            }
            else {

                matrices = Eigen.symmetricEigenvectors(F);

                E = matrices[1].diag();

                C = matrices[0].transpose();

                densityMatrix = calculateDensityMatrix(C).mmul(1 - damp).add(olddensity.mmul(damp));

            }


            numIt++;

            //System.err.println ("DIIS Error: " + DIISError);

        }




        System.out.println(moleculeName + " SCF completed: " + numIt + " iterations used");

        System.err.println ("DIIS took: " + sw.getTime());

        double e = 0;

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = 0; k < orbitals.length; k++) {
                e += 0.5 * densityMatrix.get(j, k) * (H.get(j, k) + F.get(j, k));
            }
        }

//        double checksum = 0;
//
//        for (int a = 0; a < atoms.length; a++) {
//            checksum += E(a, index);
//        }
//
//        for (int a = 0; a < atoms.length; a++) {
//            for (int b = a + 1; b < atoms.length; b++) {
//                checksum += E(a, b, index);
//            }
//        }
//
//        if (Math.abs(checksum - e) > 1E-5 || checksum != checksum) {
//            System.err.println ("I knew it!");
//            System.err.println (checksum);
//            System.err.println (e);
//            System.exit(0);
//        }


        double heat = 0;

        for (int j = 0; j < atoms.length; j++) {
            heat += atoms[j].getHeat() - atoms[j].getEisol();
            for (int k = j + 1; k < atoms.length; k++) {
                e += atoms[j].crf(atoms[k]);
            }
        }

        energy = e;

        heat += e;

        this.hf = heat / 4.3363E-2;

        if (nElectrons > 0) {
            this.homo = E.get(nElectrons / 2 - 1, 0);
        } else {
            this.homo = 0;
        }

        this.lumo = E.get(nElectrons / 2, 0);

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

        for (NDDOAtom atom : atoms) {
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


        hybridip = new double[]{0, 0, 0};

        for (int j = 0; j < atoms.length; j++) {

            if (index[j][1] != -1) {//exclude hydrogen
                hybridip[0] = hybridip[0] - 2.5416 * 2 * atoms[j].D1 * densityMatrix.get(index[j][0], index[j][1]);
                hybridip[1] = hybridip[1] - 2.5416 * 2 * atoms[j].D1 * densityMatrix.get(index[j][0], index[j][2]);
                hybridip[2] = hybridip[2] - 2.5416 * 2 * atoms[j].D1 * densityMatrix.get(index[j][0], index[j][3]);
            }
        }


        dipoletot = new double[]{chargedip[0] + hybridip[0], chargedip[1] + hybridip[1], chargedip[2] + hybridip[2]};


        dipole = Math.sqrt(dipoletot[0] * dipoletot[0] + dipoletot[1] * dipoletot[1] + dipoletot[2] * dipoletot[2]);


    }

    private double E (int atomnum, int[][] index) {

        double e = 0;

        for (int i: index[atomnum]) {
            if (i > -1) {
                e += densityMatrix.get (i, i) * orbitals[i].U();
            }
        }

        for (int i: index[atomnum]) {
            for (int j: index[atomnum]) {
                for (int k: index[atomnum]) {
                    for (int l: index[atomnum]) {
                        if (i != -1 && j != -1 && k != -1 && l != -1) {
                            e += 0.5 * densityMatrix.get(i, j) *
                                    (densityMatrix.get(k, l) * NDDO6G.OneCenterERI(orbitals[i], orbitals[j], orbitals[k], orbitals[l])
                                            - 0.5 * densityMatrix.get(k, l) * NDDO6G.OneCenterERI(orbitals[i], orbitals[k], orbitals[j], orbitals[l]));
                        }
                    }
                }
            }
        }

        return e;
    }

    private double E( int atomnum1, int atomnum2, int[][] index) {

        double e = 0;

        for (int i: index[atomnum1]) {
            for (int j: index[atomnum1]) {
                if (i != -1 && j != -1) {
                    e += densityMatrix.get(i, j) * atoms[atomnum2].V(orbitals[i], orbitals[j]);
                }
            }
        }

        for (int k: index[atomnum2]) {
            for (int l: index[atomnum2]) {
                if (k != -1 && l != -1) {
                    e += densityMatrix.get(k, l) * atoms[atomnum1].V(orbitals[k], orbitals[l]);
                }
            }
        }

        for (int i: index[atomnum1]) {
            for (int k: index[atomnum2]) {
                if (i != -1 && k != -1) {
                    e += 2 * densityMatrix.get(i, k) * NDDO6G.beta(orbitals[i], orbitals[k]);
                }
            }
        }

        for (int i: index[atomnum1]) {
            for (int j: index[atomnum1]) {
                for (int k: index[atomnum2]) {
                    for (int l: index[atomnum2]) {
                        if (i != -1 && j != -1 && k != -1 && l != -1) {
                            e +=(densityMatrix.get(i, j) * densityMatrix.get(k, l) - densityMatrix.get (i, k) * 0.5 * densityMatrix.get(j, l))
                                    * NDDO6G.getG(orbitals[i], orbitals[j], orbitals[k], orbitals[l]);
                        }
                    }
                }
            }
        }

        return e;


    }

    public DoubleMatrix getE() {
        return E;
    }


    private DoubleMatrix calculateDensityMatrix(DoubleMatrix c) {//density matrix construction by definition.
        DoubleMatrix densityMatrix = DoubleMatrix.zeros(orbitals.length, orbitals.length);
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
                densityMatrix.put(i, j, sum);
            }
        }
        return densityMatrix;
    }

    private static DoubleMatrix commutator (DoubleMatrix F, DoubleMatrix D) {

        return F.mmul(D).sub(D.mmul(F));
    }


    @Override
    public DoubleMatrix alphaDensity() {
        return this.densityMatrix.mmul(0.5);
    }

    @Override
    public DoubleMatrix betaDensity() {
        return this.densityMatrix.mmul(0.5);
    }

    @Override
    public DoubleMatrix densityMatrix() {
        return this.densityMatrix;
    }
}