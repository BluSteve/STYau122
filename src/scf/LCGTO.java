package scf;

import java.util.Arrays;

public class LCGTO {
    //LCGTO = Linear Combinations (of) GTOs i.e. contracted basis functions
    //Each LCGTO is described by an array of GTOs and an array of contraction coefficients

    protected double N;//normalisation coeffecient, should be very close to 1.
    protected double[] coordinates;
    protected int i, j, k;//see definition under GTO class
    protected GTO[] GaussArray;
    protected int n;//number of Gaussians in this contraction
    protected double[] CoeffecientArray;
    protected int shell;//principal quantum number
    protected int Z;//atomic charge
    protected int L;//angular momentum (L = 0 for s orbital, L = 1 for p orbital etc.)
    protected double zeta;
    private double[] gaussexponents;
    protected double[] coeffarray;
    protected AtomFixed atom;

    public LCGTO(double[] e, double[] c, int i, int j, int k, double[] coords, int shell, int Z, AtomFixed atom) {// e = exponent array, c = coeffecient array
        this.coeffarray = c;
        this.setGaussExponents(e);
        this.n = c.length;
        this.CoeffecientArray = c;
        this.i = i;
        this.j = j;
        this.k = k;
        this.L = i + j + k;
        this.shell = shell;
        this.coordinates = coords;
        this.Z = Z;
        this.atom = atom;

        GaussArray = new GTO[n];

        for (int a = 0; a < n; a++) {

            GTO g = new GTO(this.i, this.j, this.k, e[a], this.coordinates);
            GaussArray[a] = g;
        }


        double sum = 0;

        for (int a = 0; a < n; a++) {
            for (int b = 0; b < n; b++) {
                sum += c[a] * c[b] * GaussArray[a].getN() * GaussArray[b].getN() / Math.pow(e[a] + e[b], L + 1.5);
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
        return this.GaussArray;
    }

    public int getn() {
        return n;
    }

    public double[] getCoeffArray() {
        return CoeffecientArray;
    }

    public double[] getCoords() {
        return this.coordinates;
    }


    public boolean equals(Object b) {

        if (b instanceof LCGTO) {
            LCGTO a = (LCGTO) b;
            if (this.N == a.N && Arrays.equals(a.coordinates, coordinates) && a.i == i && a.j == j && a.k == k && Arrays.equals(a.CoeffecientArray, CoeffecientArray) && Arrays.deepEquals(a.GaussArray, GaussArray) && shell == a.shell) {

                return true;
            }

        }

        return false;
    }

    public double getZeta() {
        return this.zeta;
    }

    public boolean sameAtom(LCGTO b) {

        if (this.Z == b.Z && Arrays.equals(this.coordinates, b.coordinates)) {
            return true;
        } else {
            return false;
        }
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

    public int getshell() {

        return this.shell;
    }

    public double getN() {
        return this.N;
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


    public double[] getGaussExponents() {
        return gaussexponents;
    }


    public void setGaussExponents(double[] gaussexponents) {
        this.gaussexponents = gaussexponents;
    }

    public static double getSderiv(LCGTO X1, LCGTO X2, int tau) {
        double S = 0;

        for (int i = 0; i < X1.getn(); i++) {
            for (int j = 0; j < X2.getn(); j++) {

                S += X1.getCoeffArray()[i] * X2.getCoeffArray()[j] *
                        GTO.getSderiv(X1.getGaussArray()[i], X2.getGaussArray()[j], tau);
            }

        }

        return S * X1.getN() * X2.getN();
    }
}
