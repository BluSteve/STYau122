package scf;

import java.util.Arrays;

//this is the GTO class, and all the relevant GTO integral routines are implemented here.
//The GTO object represents a GTO of the form (x-coordinates[0])^i(y-coordinates[1])^j(z-coordinates[2])^k exp(-exponent * r^2)
public class GTO {

    private int i, j, k;// angular terms
    private double N;
    private int L; // L = i + j + k
    private double exponent; //radial exponent
    private double[] coordinates; //coordinates

    public GTO(int i, int j, int k, double exponent, double[] coordinates) {
        this.i = i;
        this.j = j;
        this.k = k;
        this.L = i + j + k;
        this.exponent = exponent;
        this.coordinates = coordinates;

        this.N = Math.pow(2 / Math.PI, 0.75) * Math.pow(2, L) * Math.pow(exponent, (2 * L + 3) / 4.0) / Math.sqrt(fact2(2 * i - 1) * fact2(2 * j - 1) * fact2(2 * k - 1));
    }

    static int fact2(int num) {
        int sum = 1;
        for (int i = num; i > 0; i -= 2) {
            sum *= i;
        }

        return sum;
    }


    public int geti() {
        return i;
    }

    public int getj() {
        return j;
    }

    public int getk() {
        return k;
    }

    public double getexponent() {
        return exponent;
    }

    public double getX() {
        return coordinates[0];
    }

    public double getY() {
        return coordinates[1];
    }

    public double getZ() {
        return coordinates[2];
    }

    public double[] getcoords() {
        return this.coordinates;
    }

    public double getN() {
        return this.N;
    }

    public int getL() {
        return this.L;
    }

    public boolean equals(Object b) {

        if (b instanceof GTO) {

            GTO c = (GTO) b;
            if (this.i == c.i && this.j == c.j && this.k == c.k && this.exponent == c.exponent && Arrays.equals(this.coordinates, c.coordinates)) {
                return true;
            }
        }
        return false;
    }


    private static double getE(int i, int j, int t, double alpha, double beta, double coordA, double coordB) { // auxillary overlap integral; recursive definition.

        double Q = coordA - coordB;
        double p = beta + alpha;
        double q = (alpha * beta) / p;
        double K = Math.pow(Math.E, (-q * Q * Q));

        if (t < 0 | t > (i + j)) {//out of bound case
            return 0;
        } else if (i == j && j == t && i == 0) {// definition
            return K;
        } else if (j == 0) {
            return ((1 / (2 * p)) * getE(i - 1, j, t - 1, alpha, beta, coordA, coordB) - (q * Q / alpha) * getE(i - 1, j, t, alpha, beta, coordA, coordB) + (t + 1) * getE(i - 1, j, t + 1, alpha, beta, coordA, coordB));
        } else {
            return ((1 / (2 * p)) * getE(i, j - 1, t - 1, alpha, beta, coordA, coordB) + (q * Q / beta) * getE(i, j - 1, t, alpha, beta, coordA, coordB) + (t + 1) * getE(i, j - 1, t + 1, alpha, beta, coordA, coordB));
        }
    }

    public static double getOverlap(GTO a, GTO b) {//evaluate <a|b>, unnormalised.

        double x = getE(a.i, b.i, 0, a.exponent, b.exponent, a.getX(), b.getX());
        double y = getE(a.j, b.j, 0, a.exponent, b.exponent, a.getY(), b.getY());
        double z = getE(a.k, b.k, 0, a.exponent, b.exponent, a.getZ(), b.getZ());


        return x * y * z * Math.pow(Math.PI / (a.exponent + b.exponent), 1.5);
    }

    public static double R(double[] P, double[] C) {

        double val = (P[0] - C[0]) * (P[0] - C[0]) + (P[1] - C[1]) * (P[1] - C[1]) + (P[2] - C[2]) * (P[2] - C[2]);

        return Math.sqrt(val);


    }

    public static double getT(GTO a, GTO b) {// evaluate <a|T|b>, unnormalised.

        double X = b.exponent * (2 * b.L + 3) * getOverlap(a, b);

        double Y = -2 * b.exponent * b.exponent *
                (getOverlap(a, new GTO(b.i + 2, b.j, b.k, b.exponent, b.coordinates)) +
                        getOverlap(a, new GTO(b.i, b.j + 2, b.k, b.exponent, b.coordinates)) +
                        getOverlap(a, new GTO(b.geti(), b.getj(), b.getk() + 2, b.getexponent(), b.getcoords())));

        double Z = -0.5 * (
                b.geti() * (b.geti() - 1) * getOverlap(a, new GTO(b.geti() - 2, b.getj(), b.getk(), b.getexponent(), b.getcoords())) +
                        b.getj() * (b.getj() - 1) * getOverlap(a, new GTO(b.geti(), b.getj() - 2, b.getk(), b.getexponent(), b.getcoords())) +
                        b.getk() * (b.getk() - 1) * getOverlap(a, new GTO(b.geti(), b.getj(), b.getk() - 2, b.getexponent(), b.getcoords())));

        return (X + Y + Z);
    }

    private static double getR(int t, int u, int v, int n, double p, double[] P, double R) {// auxillary potential energy integral; recursive definition.

        double T = p * R * R;
        double returnVal = 0;

        if (t == 0 && u == 0 && v == 0) {
            returnVal += Math.pow(-2 * p, n) * Boys(n, T);//base case
        } else if (t == 0 && u == 0) {

            if (v > 1) {
                returnVal += (v - 1) * getR(t, u, v - 2, n + 1, p, P, R);
            }
            returnVal += P[2] * getR(t, u, v - 1, n + 1, p, P, R);
        } else if (t == 0) {

            if (u > 1) {
                returnVal += (u - 1) * getR(t, u - 2, v, n + 1, p, P, R);
            }
            returnVal += P[1] * getR(t, u - 1, v, n + 1, p, P, R);

        } else {

            if (t > 1) {
                returnVal += (t - 1) * getR(t - 2, u, v, n + 1, p, P, R);
            }
            returnVal += P[0] * getR(t - 1, u, v, n + 1, p, P, R);

        }

        return returnVal;
    }

    private static double Boys(int n, double T) {//Boys function Fn(T). Used in the potential energy and two-electron integrals.
        return SpecialFunctions.confluentHypergeometric1F1(n + 0.5, n + 1.5, -T) / (2 * n + 1);
    }

    public static double erf2(double z) {
        double t = 1.0 / (1.0 + 0.47047 * Math.abs(z));
        double poly = t * (0.3480242 + t * (-0.0958798 + t * (0.7478556)));
        double ans = 1.0 - poly * Math.exp(-z * z);
        if (z >= 0) return ans;
        else return -ans;
    }

    public static double getV(GTO a, GTO b, double[] C) {//evaluate <a|V|b>. This is not normalized.

        double V = 0;
        double[] P = new double[3];

        for (int i = 0; i < 3; i++) {
            P[i] = (a.getexponent() * a.getcoords()[i] + b.getexponent() * b.getcoords()[i]) / (a.getexponent() + b.getexponent());
        }

        double R = R(P, C);
        double p = a.getexponent() + b.getexponent();

        for (int i = 0; i < (a.geti() + b.geti() + 1); i++) {
            for (int j = 0; j < a.getj() + b.getj() + 1; j++) {
                for (int k = 0; k < a.getk() + b.getk() + 1; k++) {

                    V += getE(a.geti(), b.geti(), i, a.getexponent(), b.getexponent(), a.getX(), b.getX()) *
                            getE(a.getj(), b.getj(), j, a.getexponent(), b.getexponent(), a.getY(), b.getY()) *
                            getE(a.getk(), b.getk(), k, a.getexponent(), b.getexponent(), a.getZ(), b.getZ()) *
                            getR(i, j, k, 0, p, new double[]{P[0] - C[0], P[1] - C[1], P[2] - C[2]}, R);
                }
            }
        }

        return V * 2 * Math.PI / p;
    }

    public static double getG(GTO a, GTO b, GTO c, GTO d) {//evaluate (ab|cd) (Chemist notation). This is not normalized.

        double integral = 0;
        double[] P = new double[3];

        for (int i = 0; i < 3; i++) {
            P[i] = (a.getexponent() * a.getcoords()[i] + b.getexponent() * b.getcoords()[i]) / (a.getexponent() + b.getexponent());
        }

        double[] C = new double[3];

        for (int i = 0; i < 3; i++) {

            C[i] = (c.getexponent() * c.getcoords()[i] + d.getexponent() * d.getcoords()[i]) / (c.getexponent() + d.getexponent());
        }

        double[] Q = new double[]{P[0] - C[0], P[1] - C[1], P[2] - C[2]};

        double p = ((a.getexponent() + b.getexponent()) * (c.getexponent() + d.getexponent())) / (a.getexponent() + b.getexponent() + c.getexponent() + d.getexponent());

        for (int i = 0; i < (a.geti() + b.geti() + 1); i++) {
            for (int j = 0; j < (a.getj() + b.getj() + 1); j++) {
                for (int k = 0; k < (a.getk() + b.getk() + 1); k++) {
                    for (int nu = 0; nu < (c.geti() + d.geti() + 1); nu++) {
                        for (int kappa = 0; kappa < (c.getj() + d.getj() + 1); kappa++) {
                            for (int phi = 0; phi < (c.getk() + d.getk() + 1); phi++) {

                                integral += getE(a.geti(), b.geti(), i, a.getexponent(), b.getexponent(), a.getX(), b.getX()) *
                                        getE(a.getj(), b.getj(), j, a.getexponent(), b.getexponent(), a.getY(), b.getY()) *
                                        getE(a.getk(), b.getk(), k, a.getexponent(), b.getexponent(), a.getZ(), b.getZ()) *
                                        getE(c.geti(), d.geti(), nu, c.getexponent(), d.getexponent(), c.getX(), d.getX()) *
                                        getE(c.getj(), d.getj(), kappa, c.getexponent(), d.getexponent(), c.getY(), d.getY()) *
                                        getE(c.getk(), d.getk(), phi, c.getexponent(), d.getexponent(), c.getZ(), d.getZ()) *
                                        Math.pow(-1, nu + kappa + phi) *
                                        getR(i + nu, j + kappa, k + phi, 0, p, Q, R(P, C));

                            }
                        }
                    }
                }
            }
        }

        integral *= (2 * Math.pow(Math.PI, 2.5) / ((a.getexponent() + b.getexponent()) * (c.getexponent() + d.getexponent()) *
                Math.sqrt(a.getexponent() + b.getexponent() + c.getexponent() + d.getexponent())));

        return integral;
    }

    public static double getS(GTO a, GTO b) {

        double xdir = 0, ydir = 0, zdir = 0;

        double[] Q = new double[]{a.coordinates[0] - b.coordinates[0], a.coordinates[1] - b.coordinates[1], a.coordinates[2] - b.coordinates[2]};
        double p = a.exponent + b.exponent;
        double q = (a.exponent * b.exponent) / p;
        double[] k = new double[]{Math.exp(-q * Q[0] * Q[0]), Math.exp(-q * Q[1] * Q[1]), Math.exp(-q * Q[2] * Q[2])};

        if (a.i + b.i == 0) {
            xdir = k[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.i == 0 && b.i == 1) {
            xdir = a.exponent / (a.exponent + b.exponent) * Q[0] * k[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.i == 1 && b.i == 0) {
            xdir = -b.exponent / (a.exponent + b.exponent) * Q[0] * k[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.i == 1 && b.i == 1) {
            xdir = k[0] * Math.sqrt(Math.PI) / (2 * Math.pow(a.exponent + b.exponent, 1.5)) - k[0] * q / p * Q[0] * Q[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        }

        if (a.j + b.j == 0) {
            ydir = k[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.j == 0 && b.j == 1) {
            ydir = a.exponent / (a.exponent + b.exponent) * Q[1] * k[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.j == 1 && b.j == 0) {
            ydir = -b.exponent / (a.exponent + b.exponent) * Q[1] * k[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.j == 1 && b.j == 1) {
            ydir = k[1] * Math.sqrt(Math.PI) / (2 * Math.pow(a.exponent + b.exponent, 1.5)) - k[1] * q / p * Q[1] * Q[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        }

        if (a.k + b.k == 0) {
            zdir = k[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.k == 0 && b.k == 1) {
            zdir = a.exponent / (a.exponent + b.exponent) * Q[2] * k[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.k == 1 && b.k == 0) {
            zdir = -b.exponent / (a.exponent + b.exponent) * Q[2] * k[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.k == 1 && b.k == 1) {
            zdir = k[2] * Math.sqrt(Math.PI) / (2 * Math.pow(a.exponent + b.exponent, 1.5)) - k[2] * q / p * Q[2] * Q[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        }

        return xdir * ydir * zdir * a.N * b.N;
    }

    public static double getSderiv(GTO a, GTO b, int tau) {
        double[] Q = new double[]{a.coordinates[0] - b.coordinates[0], a.coordinates[1] - b.coordinates[1], a.coordinates[2] - b.coordinates[2]};
        double p = a.exponent + b.exponent;
        double q = (a.exponent * b.exponent) / p;
        double[] k = new double[]{Math.exp(-q * Q[0] * Q[0]), Math.exp(-q * Q[1] * Q[1]), Math.exp(-q * Q[2] * Q[2])};
        double[] kprime = new double[]{-2 * q * Q[0] * k[0], -2 * q * Q[1] * k[1], -2 * q * Q[2] * k[2]};
        double[] Qprime = new double[]{1.0, 1.0, 1.0};

        double xderiv = 0, yderiv = 0, zderiv = 0;

        double xdir = 0, ydir = 0, zdir = 0;

        if (a.i + b.i == 0) {
            xdir = k[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            xderiv = kprime[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);

        } else if (a.i == 0 && b.i == 1) {
            xdir = a.exponent / (a.exponent + b.exponent) * Q[0] * k[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            xderiv = a.exponent / (a.exponent + b.exponent) * Qprime[0] * k[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent)
                    + a.exponent / (a.exponent + b.exponent) * Q[0] * kprime[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.i == 1 && b.i == 0) {
            xdir = -b.exponent / (a.exponent + b.exponent) * Q[0] * k[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            xderiv = -b.exponent / (a.exponent + b.exponent) * Qprime[0] * k[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent)
                    - b.exponent / (a.exponent + b.exponent) * Q[0] * kprime[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.i == 1 && b.i == 1) {
            xdir = k[0] * Math.sqrt(Math.PI) / (2 * Math.pow(a.exponent + b.exponent, 1.5)) - k[0] * q / p * Q[0] * Q[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            xderiv = kprime[0] * Math.sqrt(Math.PI) / (2 * Math.pow(a.exponent + b.exponent, 1.5))
                    - kprime[0] * q / p * Q[0] * Q[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent)
                    - k[0] * q / p * 2 * Qprime[0] * Q[0] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        }

        if (a.j + b.j == 0) {
            ydir = k[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            yderiv = kprime[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.j == 0 && b.j == 1) {
            ydir = a.exponent / (a.exponent + b.exponent) * Q[1] * k[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            yderiv = a.exponent / (a.exponent + b.exponent) * Qprime[1] * k[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent)
                    + a.exponent / (a.exponent + b.exponent) * Q[1] * kprime[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.j == 1 && b.j == 0) {
            ydir = -b.exponent / (a.exponent + b.exponent) * Q[1] * k[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            yderiv = -b.exponent / (a.exponent + b.exponent) * Qprime[1] * k[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent)
                    - b.exponent / (a.exponent + b.exponent) * Q[1] * kprime[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);

        } else if (a.j == 1 && b.j == 1) {
            ydir = k[1] * Math.sqrt(Math.PI) / (2 * Math.pow(a.exponent + b.exponent, 1.5)) - k[1] * q / p * Q[1] * Q[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            yderiv = kprime[1] * Math.sqrt(Math.PI) / (2 * Math.pow(a.exponent + b.exponent, 1.5))
                    - kprime[1] * q / p * Q[1] * Q[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent)
                    - k[1] * q / p * 2 * Qprime[1] * Q[1] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);

        }

        if (a.k + b.k == 0) {
            zdir = k[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            zderiv = kprime[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.k == 0 && b.k == 1) {
            zdir = a.exponent / (a.exponent + b.exponent) * Q[2] * k[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            zderiv = a.exponent / (a.exponent + b.exponent) * Qprime[2] * k[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent)
                    + a.exponent / (a.exponent + b.exponent) * Q[2] * kprime[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.k == 1 && b.k == 0) {
            zdir = -b.exponent / (a.exponent + b.exponent) * Q[2] * k[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            zderiv = -b.exponent / (a.exponent + b.exponent) * Qprime[2] * k[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent)
                    - b.exponent / (a.exponent + b.exponent) * Q[2] * kprime[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        } else if (a.k == 1 && b.k == 1) {
            zdir = k[2] * Math.sqrt(Math.PI) / (2 * Math.pow(a.exponent + b.exponent, 1.5)) - k[2] * q / p * Q[2] * Q[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
            zderiv = kprime[2] * Math.sqrt(Math.PI) / (2 * Math.pow(a.exponent + b.exponent, 1.5))
                    - kprime[2] * q / p * Q[2] * Q[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent)
                    - k[2] * q / p * 2 * Qprime[2] * Q[2] * Math.sqrt(Math.PI) / Math.sqrt(a.exponent + b.exponent);
        }

        switch (tau) {
            case 0:
                return xderiv * ydir * zdir * a.N * b.N;
            case 1:
                return xdir * yderiv * zdir * a.N * b.N;
            case 2:
                return xdir * ydir * zderiv * a.N * b.N;

        }
        return 0;
    }

    public static double getSderivfinite(GTO a, GTO b, int tau) {

        double[] newcoords = a.coordinates.clone();
        newcoords[tau] += 1E-8;

        GTO newa = new GTO(a.i, a.j, a.k, a.exponent, newcoords);

        return 1E8 * (GTO.getS(newa, b) - GTO.getS(a, b));

    }
}
