package scf;

import java.util.Arrays;

public class LCGTO { // HAS AtomFixed, OrbitalProperties, c, e
    //LCGTO = Linear Combinations (of) GTOs i.e. contracted basis functions
    //Each LCGTO is described by an array of GTOs and an array of contraction coefficients

    protected double N;//normalisation coefficient, should be very close to 1.
    protected double[] coordinates;
    protected int i, j, k;//see definition under GTO class
    protected GTO[] gaussArray;
    protected int n;//number of Gaussians in this contraction
    protected double[] coefficientArray;
    protected int shell;//principal quantum number
    protected int Z;//atomic charge
    protected int L;//angular momentum (L = 0 for s orbital, L = 1 for p orbital etc.)
    protected double zeta;
    private double[] gaussExponents;
    protected AtomFixed atom;

    public LCGTO(double[] e, double[] c, AtomFixed atom, OrbitalProperties orbital) {// e = exponent array, c = coefficient array
        this.gaussExponents = e;
        this.coefficientArray = c;
        this.n = c.length;
        this.i = orbital.getConfig()[0];
        this.j = orbital.getConfig()[1];
        this.k = orbital.getConfig()[2];
        this.L = i + j + k;
        this.shell = orbital.getShell();
        this.coordinates = atom.getCoordinates();
        this.Z = atom.getAtomProperties().getZ();
        this.atom = atom;

        gaussArray = new GTO[n];

        for (int a = 0; a < n; a++) {
            GTO g = new GTO(this.i, this.j, this.k, e[a], this.coordinates);
            gaussArray[a] = g;
        }

        double sum = 0;

        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                sum += c[a] * c[b] * gaussArray[a].getN() * gaussArray[b].getN() / Math.pow(e[a] + e[b], L + 1.5);
            }
        }

        sum *= Math.pow(Math.PI, 1.5) * GTO.fact2(2 * i - 1) * GTO.fact2(2 * j - 1) * GTO.fact2(2 * k - 1) / Math.pow(2, L);
        sum = Math.pow(sum, -0.5);
        N = sum;
    }

    public LCGTO() {

    }

    public int getZ() {
        return this.Z;
    }

    public AtomFixed getAtom() {
        return this.atom;
    }

    public GTO[] getGaussArray() {
        return this.gaussArray;
    }

    public int getn() {
        return n;
    }

    public double[] getCoeffArray() {
        return coefficientArray;
    }

    public double[] getCoords() {
        return this.coordinates;
    }

    public double getZeta() {
        return this.zeta;
    }

    public int getL() {
        return i + j + k;
    }

    public int geti() {
        return i;
    }

    public int getj() {
        return this.j;
    }

    public int getk() {
        return this.k;
    }

    public int getShell() {
        return this.shell;
    }

    public double getN() {
        return this.N;
    }

    public double[] getGaussExponents() {
        return gaussExponents;
    }

    public static double getS(LCGTO X1, LCGTO X2) {//normalised overlap integral
        double S = 0;

        for (int i = 0; i < X1.getn(); i++) {
            for (int j = 0; j < X2.getn(); j++) {
                S += X1.getGaussArray()[i].getN() * X2.getGaussArray()[j].getN() *
                        X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
                        GTO.getOverlap(X1.getGaussArray()[i], X2.getGaussArray()[j]);
            }

        }

        return S * X1.getN() * X2.getN();
    }

    public static double getT(LCGTO X1, LCGTO X2) {//normalised kinetic energy integral
        double T = 0;

        for (int i = 0; i < X1.getn(); i++) {
            for (int j = 0; j < X2.getn(); j++) {
                T += X1.getGaussArray()[i].getN() * X2.getGaussArray()[j].getN() *
                        X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
                        GTO.getT(X1.getGaussArray()[i], X2.getGaussArray()[j]);
            }
        }

        return T * X1.getN() * X2.getN();

    }

    public static double getV(LCGTO X1, LCGTO X2, double[] coords) {//normalised nuclear-electron potential energy integral
        double V = 0;

        for (int i = 0; i < X1.getn(); i++) {
            for (int j = 0; j < X2.getn(); j++) {
                V += X1.getGaussArray()[i].getN() * X2.getGaussArray()[j].getN() *
                        X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
                        GTO.getV(X1.getGaussArray()[i], X2.getGaussArray()[j], coords);
            }
        }

        return V * X1.getN() * X2.getN();
    }

    public static double getG(LCGTO X1, LCGTO X2, LCGTO X3, LCGTO X4) {//normalised electron-repulsion integral
        double integral = 0;

        for (int i = 0; i < X1.getn(); i++) {
            for (int j = 0; j < X2.getn(); j++) {
                for (int k = 0; k < X3.getn(); k++) {
                    for (int l = 0; l < X4.getn(); l++) {
                        integral += X1.getGaussArray()[i].getN() * X2.getGaussArray()[j].getN() * X3.getGaussArray()[k].getN() * X4.getGaussArray()[l].getN() *
                                X1.getCoeffArray()[i] * X2.getCoeffArray()[j] * X3.getCoeffArray()[k] * X4.getCoeffArray()[l] *
                                GTO.getG(X1.getGaussArray()[i], X2.getGaussArray()[j], X3.getGaussArray()[k], X4.getGaussArray()[l]);
                    }
                }
            }
        }

        return integral * X1.getN() * X2.getN() * X3.getN() * X4.getN();
    }

    public static double getSDeriv(LCGTO X1, LCGTO X2, int tau) {
        double S = 0;

        for (int i = 0; i < X1.getn(); i++) {
            for (int j = 0; j < X2.getn(); j++) {
                S += X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
                        GTO.getSderiv(X1.getGaussArray()[i], X2.getGaussArray()[j], tau);
            }

        }

        return S * X1.getN() * X2.getN();
    }

    public boolean equals(Object b) {
        if (b instanceof LCGTO) {
            LCGTO a = (LCGTO) b;
            if (this.N == a.N && Arrays.equals(a.coordinates, coordinates) && a.i == i && a.j == j && a.k == k && Arrays.equals(a.coefficientArray, coefficientArray) && Arrays.deepEquals(a.gaussArray, gaussArray) && shell == a.shell) {
                return true;
            }

        }
        return false;
    }
}
