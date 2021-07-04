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
    protected AtomFixed atom;
    private double[] gaussExponents;

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

    public static double getS(LCGTO X1, LCGTO X2) {//normalised overlap integral
        double S = 0;

        for (int i = 0; i < X1.getn(); i++) {
            for (int j = 0; j < X2.getn(); j++) {
                S += X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
                        GTO.getS(X1.getGaussArray()[i], X2.getGaussArray()[j]);
            }

        }

        return S * X1.getN() * X2.getN();
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

    public static double getSDeriv2(LCGTO X1, LCGTO X2, int tau1, int tau2) {
        double S = 0;

        for (int i = 0; i < X1.getn(); i++) {
            for (int j = 0; j < X2.getn(); j++) {
                S += X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
                        GTO.getSderiv2(X1.getGaussArray()[i], X2.getGaussArray()[j], tau1, tau2);
            }

        }

        return S * X1.getN() * X2.getN();
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

    public boolean equals(Object b) {
        if (b instanceof LCGTO) {
            LCGTO a = (LCGTO) b;
            return this.N == a.N && Arrays.equals(a.coordinates, coordinates) && a.i == i && a.j == j && a.k == k && Arrays.equals(a.coefficientArray, coefficientArray) && Arrays.deepEquals(a.gaussArray, gaussArray) && shell == a.shell;

        }
        return false;
    }
}
