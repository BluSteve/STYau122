package mndoparam.mndo;

import org.jblas.DoubleMatrix;

import java.util.Arrays;

public abstract class MNDOSolution {
    public double energy, homo, lumo, hf, dipole;
    public double[] chargedip, hybridip, dipoletot;
    public int charge, multiplicity;

    protected int[][] missingIndex, index;
    protected int[] atomNumber;
    protected double damp = 0.8;
    protected int nElectrons;
    protected DoubleMatrix H;
    protected MNDO6G[] orbitals;

    public MNDOSolution(MNDOAtom[] atoms, int charge) {
        for (MNDOAtom a : atoms) {
            nElectrons += a.getAtomProperties().getQ();
        }
        nElectrons -= charge;

        this.charge = charge;
        int i = 0;

        for (MNDOAtom a : atoms) {
            i += a.getOrbitals().length;
        }

        orbitals = new MNDO6G[i];

        i = 0;

        index = new int[atoms.length][4];
        atomNumber = new int[orbitals.length];
        int count = 0;
        int count2;
        for (MNDOAtom a : atoms) {
            count2 = 0;
            for (MNDO6G orbital : a.getOrbitals()) {
                orbitals[i] = orbital;
                index[count][count2] = i;
                atomNumber[i] = count;
                i++;
                count2++;
            }


            if (a.getAtomProperties().getZ() == 1) {
                index[count][1] = -1;
                index[count][2] = -1;
                index[count][3] = -1;

            }
            count++;
        }

        missingIndex = new int[atoms.length][4 * atoms.length - 4];

        for (int j = 0; j < atoms.length; j++) {
            for (int k = 0; k < 4 * atoms.length - 4; k++) {
                missingIndex[j][k] = -1;
            }
        }

        for (int j = 0; j < atoms.length; j++) {
            int[] nums = new int[]{index[j][0], index[j][1], index[j][2], index[j][3]};
            int counter = 0;
            for (int k = 0; k < orbitals.length; k++) {
                if (nums[0] != k && nums[1] != k && nums[2] != k && nums[3] != k) {
                    missingIndex[j][counter] = k;
                    counter++;
                }
            }
        }

        H = new DoubleMatrix(orbitals.length, orbitals.length);

        //filling up the core matrix in accordance with MNDO formalism

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                if (k == j) {
                    double Huu = orbitals[j].U();
                    for (MNDOAtom a : atoms) {
                        if (!Arrays.equals(a.getCoordinates(), orbitals[j].getCoords())) { // case 1
                            Huu += a.V(orbitals[j], orbitals[k]);
                        }
                    }
                    H.put(j, k, Huu);
                } else if (atomNumber[j] == atomNumber[k]) { // case 2
                    double Huv = 0;
                    for (MNDOAtom a : atoms) {
                        if (!Arrays.equals(a.getCoordinates(), orbitals[j].getCoords())) { // TODO remove duplicate code
                            Huv += a.V(orbitals[j], orbitals[k]);
                        }
                    }
                    H.put(j, k, Huv);
                    H.put(k, j, Huv);
                } else { // case 3
                    double Huk = MNDO6G.beta(orbitals[j], orbitals[k]);
                    H.put(j, k, Huk);
                    H.put(k, j, Huk);
                }
            }
        }
    }

    protected static boolean isSimilar(DoubleMatrix x, DoubleMatrix y, double limit) {
        for (int i = 0; i < y.rows; i++) {
            for (int j = 0; j < y.columns; j++) {
                if (Math.abs(x.get(i, j) - y.get(i, j)) > limit) {
                    return false;
                }
            }
        }
        return true;
    }

    public abstract DoubleMatrix alphaDensity();

    public abstract DoubleMatrix betaDensity();

    public abstract DoubleMatrix densityMatrix();

}
