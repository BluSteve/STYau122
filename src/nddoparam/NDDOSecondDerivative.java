package nddoparam;

import org.jblas.DoubleMatrix;
import org.jblas.Singular;

import java.util.*;

import scf.GTO;
import scf.Utils;


public class NDDOSecondDerivative {


    private static double generalizedform(double a, double b, double R) {
        double denom = (R + a) * (R + a) + b;
        return 1 / (R * R * Math.pow(denom, 1.5)) * (a / R + 3 * (R + a) * (R + a) / (denom));
    }

    private static double generalizedform2(double a, double b, double R) {
        return -(R + a) / (Math.pow((R + a) * (R + a) + b, 1.5) * R);
    }

    private static double qqderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a00 = p01 + p02;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {

            sum = generalizedform2(0, a00 * a00, R);
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a00 * a00, R);
    }

    private static double quzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a01 = p01 + p12;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {
            sum = generalizedform2(D12, a01 * a01, R) * 0.5 + generalizedform2(-D12, a01 * a01, R) * -0.5;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(D12, a01 * a01, R) * 0.5
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D12, a01 * a01, R) * -0.5;
    }

    private static double qQpipideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a02 = p01 + p22;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {
            sum = generalizedform2(0, 4 * D22 * D22 + a02 * a02, R) * 0.5 + generalizedform2(0, a02 * a02, R) * -0.5;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D22 * D22 + a02 * a02, R) * 0.5
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a02 * a02, R) * -0.5;
    }

    private static double qQzzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a02 = p01 + p22;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {
            sum = generalizedform2(2 * D22, a02 * a02, R) * 0.25 + generalizedform2(-2 * D22, a02 * a02, R) * 0.25 + generalizedform2(0, a02 * a02, R) * -0.5;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(2 * D22, a02 * a02, R) * 0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D22, a02 * a02, R) * 0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a02 * a02, R) * -0.5;
    }

    private static double upiupideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a11 = p11 + p12;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {
            sum = generalizedform2(0, (D11 - D12) * (D11 - D12) + a11 * a11, R) * 0.5 + generalizedform2(0, (D11 + D12) * (D11 + D12) + a11 * a11, R) * -0.5;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, (D11 - D12) * (D11 - D12) + a11 * a11, R) * 0.5
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, (D11 + D12) * (D11 + D12) + a11 * a11, R) * -0.5;
    }

    private static double uzuzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a11 = p11 + p12;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {
            sum = generalizedform2(D11 - D12, a11 * a11, R) * 0.25 + generalizedform2(D11 + D12, a11 * a11, R) * -0.25 + generalizedform2(-D11 - D12, a11 * a11, R) * -0.25 + generalizedform2(-D11 + D12, a11 * a11, R) * 0.25;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(D11 - D12, a11 * a11, R) * 0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(D11 + D12, a11 * a11, R) * -0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11 - D12, a11 * a11, R) * -0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11 + D12, a11 * a11, R) * 0.25;
    }

    private static double upiQpizderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a12 = p11 + p22;
        double R = GTO.R(xA, xB);
        double denom = (D11 - D22) * (D11 - D22) + a12 * a12;
        double denom2 = (D11 + D22) * (D11 + D22) + a12 * a12;
        if (tau1 == tau2) {
            sum = generalizedform2(-D22, denom, R) * -0.25 + generalizedform2(-D22, denom2, R) * 0.25 + generalizedform2(+D22, denom, R) * 0.25 + generalizedform2(+D22, denom2, R) * -0.25;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D22, denom, R) * -0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D22, denom2, R) * 0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D22, denom, R) * 0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D22, denom2, R) * -0.25;
    }

    private static double uzQpipideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a12 = p11 + p22;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {
            sum = generalizedform2(+D11, 4 * D22 * D22 + a12 * a12, R) * -0.25 + generalizedform2(-D11, 4 * D22 * D22 + a12 * a12, R) * 0.25 + generalizedform2(+D11, a12 * a12, R) * 0.25 + generalizedform2(-D11, a12 * a12, R) * -0.25;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11, 4 * D22 * D22 + a12 * a12, R) * -0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11, 4 * D22 * D22 + a12 * a12, R) * 0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11, a12 * a12, R) * 0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11, a12 * a12, R) * -0.25;
    }

    private static double uzQzzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a12 = p11 + p22;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {
            sum = generalizedform2(+D11 - 2 * D22, a12 * a12, R) * -0.125 + generalizedform2(-D11 - 2 * D22, a12 * a12, R) * 0.125 + generalizedform2(+D11 + 2 * D22, a12 * a12, R) * -0.125
                    + generalizedform2(-D11 + 2 * D22, a12 * a12, R) * 0.125 + generalizedform2(+D11, a12 * a12, R) * 0.25 + generalizedform2(-D11, a12 * a12, R) * -0.25;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11 - 2 * D22, a12 * a12, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11 - 2 * D22, a12 * a12, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11 + 2 * D22, a12 * a12, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11 + 2 * D22, a12 * a12, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D11, a12 * a12, R) * 0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D11, a12 * a12, R) * -0.25;
    }

    private static double QpipiQpipideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a22 = p21 + p22;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {
            sum = generalizedform2(0, 4 * (D21 - D22) * (D21 - D22) + a22 * a22, R) * 0.125 + generalizedform2(0, 4 * (D21 + D22) * (D21 + D22) + a22 * a22, R) * 0.125
                    + generalizedform2(0, 4 * D21 * D21 + a22 * a22, R) * -0.25 + generalizedform2(0, 4 * D22 * D22 + a22 * a22, R) * -0.25 + generalizedform2(0, a22 * a22, R) * 0.25;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * (D21 - D22) * (D21 - D22) + a22 * a22, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * (D21 + D22) * (D21 + D22) + a22 * a22, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D21 * D21 + a22 * a22, R) * -0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D22 * D22 + a22 * a22, R) * -0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a22 * a22, R) * 0.25;
    }

    private static double QxxQyyderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a22 = p21 + p22;
        double R = GTO.R(xA, xB);

        if (tau1 == tau2) {
            sum = generalizedform2(0, 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, R) * 0.25 + generalizedform2(0, 4 * D21 * D21 + a22 * a22, R) * -0.25 + generalizedform2(0, 4 * D22 * D22 + a22 * a22, R) * -0.25 + generalizedform2(0, a22 * a22, R) * 0.25;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, R) * 0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D21 * D21 + a22 * a22, R) * -0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, 4 * D22 * D22 + a22 * a22, R) * -0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a22 * a22, R) * 0.25;
    }

    private static double QpipiQzzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a22 = p21 + p22;
        double R = GTO.R(xA, xB);
        double denom = 4 * D21 * D21 + a22 * a22;
        if (tau1 == tau2) {
            sum = generalizedform2(-2 * D22, denom, R) * 0.125 + generalizedform2(+2 * D22, denom, R) * 0.125 + generalizedform2(-2 * D22, a22 * a22, R) * -0.125 + generalizedform2(+2 * D22, a22 * a22, R) * -0.125 + generalizedform2(0, denom, R) * -0.25 + generalizedform2(0, a22 * a22, R) * 0.25;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D22, denom, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D22, denom, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D22, a22 * a22, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D22, a22 * a22, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, denom, R) * -0.25
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a22 * a22, R) * 0.25;
    }


    private static double QzzQzzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a22 = p21 + p22;
        double R = GTO.R(xA, xB);
        if (tau1 == tau2) {
            sum = generalizedform2(+2 * D21 - 2 * D22, a22 * a22, R) * 0.0625 + generalizedform2(+2 * D21 + 2 * D22, a22 * a22, R) * 0.0625 + generalizedform2(-2 * D21 - 2 * D22, a22 * a22, R) * 0.0625 + generalizedform2(-2 * D21 + 2 * D22, a22 * a22, R) * 0.0625
                    + generalizedform2(+2 * D21, a22 * a22, R) * -0.125 + generalizedform2(-2 * D21, a22 * a22, R) * -0.125 + generalizedform2(+2 * D22, a22 * a22, R) * -0.125 + generalizedform2(-2 * D22, a22 * a22, R) * -0.125 + generalizedform2(0, a22 * a22, R) * 0.25;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D21 - 2 * D22, a22 * a22, R) * 0.0625
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D21 + 2 * D22, a22 * a22, R) * 0.0625
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D21 - 2 * D22, a22 * a22, R) * 0.0625
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D21 + 2 * D22, a22 * a22, R) * 0.0625
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D21, a22 * a22, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D21, a22 * a22, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+2 * D22, a22 * a22, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-2 * D22, a22 * a22, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(0, a22 * a22, R) * 0.25;
    }

    private static double QpizQpizderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        double sum = 0;

        double a22 = p21 + p22;
        double R = GTO.R(xA, xB);
        double denom1 = (D21 - D22) * (D21 - D22) + a22 * a22;
        double denom2 = (D21 + D22) * (D21 + D22) + a22 * a22;
        if (tau1 == tau2) {
            sum = generalizedform2(+D21 - D22, denom1, R) * 0.125 + generalizedform2(+D21 - D22, denom2, R) * -0.125 + generalizedform2(+D21 + D22, denom1, R) * -0.125 + generalizedform2(+D21 + D22, denom2, R) * 0.125
                    + generalizedform2(-D21 - D22, denom1, R) * -0.125 + generalizedform2(-D21 - D22, denom2, R) * 0.125 + generalizedform2(-D21 + D22, denom1, R) * 0.125 + generalizedform2(-D21 + D22, denom2, R) * -0.125;
        }
        return sum + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D21 - D22, denom1, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D21 - D22, denom2, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D21 + D22, denom1, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(+D21 + D22, denom2, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D21 - D22, denom1, R) * -0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D21 - D22, denom2, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D21 + D22, denom1, R) * 0.125
                + (xB[tau1] - xA[tau1]) * (xB[tau2] - xA[tau2]) * generalizedform(-D21 + D22, denom2, R) * -0.125;
    }

    private static double ssssderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double ssppippideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double sspzpzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double ppippissderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
    }

    private static double pzpzssderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
    }

    private static double ppippippippideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) + QpipiQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double pxpxpypyderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) + QxxQyyderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double ppippipzpzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) + QpipiQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double pzpzppippideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) + QpipiQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
    }

    private static double pzpzpzpzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return qqderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) + qQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) + QzzQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double spzssderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return -quzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
    }

    private static double spzppippideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return -quzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) + uzQpipideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double spzpzpzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return -quzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2) + uzQzzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double ssspzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return quzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double ppippispzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return quzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) - uzQpipideriv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
    }

    private static double pzpzspzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return quzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) - uzQzzderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
    }

    private static double sppisppideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return upiupideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double spzspzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return uzuzderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double sppippipzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return upiQpizderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double ppipzsppideriv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return -upiQpizderiv2(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau1, tau2);
    }

    private static double ppipzppipzderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return QpizQpizderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2);
    }

    private static double pxpypxpyderiv2(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau1, int tau2) {
        return 0.5 * (ppippippippideriv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2) - pxpxpypyderiv2(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau1, tau2));
    }

    protected static double LocalTwoCenterERIderiv2(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau1, int tau2) {

        double[] A = a.getCoords();
        double[] C = c.getCoords();
        //(??|??)
        switch (a.getL()) {

            case 0://(s?|??)

                switch (b.getL()) {

                    case 0: //(ss|??)

                        switch (c.getL()) {

                            case 0: //(ss|s?);

                                switch (d.getL()) {

                                    case 0://(ss|ss)
                                        return ssssderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);

                                    case 1:
                                        if (d.getk() == 1) {//(ss|spz)
                                            return ssspzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                        } else {//(ss|sppi) = 0
                                            return 0;
                                        }
                                    default:
                                        System.err.println("oh no");
                                        return 0;
                                }

                            case 1: //(ss|p?)
                                if (c.getk() == 1) {//(ss|pz?)

                                    switch (d.getL()) {

                                        case 0://(ss|pzs)
                                            return ssspzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);

                                        case 1:
                                            if (d.getk() == 1) {//(ss|pzpz)
                                                return sspzpzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                            } else {//(ss|pzppi) = 0
                                                return 0;
                                            }
                                        default:
                                            return 0;
                                    }
                                } else {//(ss|ppi?)

                                    if (d.getL() == 1 && d.getk() == 0 && c.geti() == d.geti() && c.getj() == d.getj()) {//(ss|ppippi)
                                        return ssppippideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                    } else {//all others are 0
                                        return 0;
                                    }
                                }
                            default:
                                System.err.println("oh no");
                                return 0;

                        }
                    case 1: //(sp|??)

                        if (b.getk() == 1) {//(spz|??)

                            switch (c.getL()) {

                                case 0://(spz|s?)

                                    switch (d.getL()) {

                                        case 0://(spz|ss)
                                            return spzssderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);

                                        case 1:
                                            if (d.getk() == 1) {//(spz|spz)
                                                return spzspzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                            } else {
                                                return 0;
                                            }
                                    }

                                case 1:
                                    if (c.getk() == 1) {//(spz|pz?)

                                        switch (d.getL()) {

                                            case 0://(spz|pzs)
                                                return spzspzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);

                                            case 1:
                                                if (d.getk() == 1) {//(spz|pzpz)
                                                    return spzpzpzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                } else {//(spz|pzppi) = 0
                                                    return 0;
                                                }
                                        }
                                    } else {//(spz|ppi?)
                                        if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
                                            return spzppippideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                        } else {
                                            return 0;
                                        }
                                    }
                                default:
                                    System.err.println("oh no");
                                    return 0;
                            }
                        } else {//(sppi|??)

                            switch (c.getL()) {
                                case 0://(sppi|s?)
                                    if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(sppi|sppi)
                                        return sppisppideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                    } else {
                                        return 0;
                                    }
                                case 1:
                                    if (c.getk() == 1) {
                                        if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(sppi|pzppi)
                                            return sppippipzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                        } else {
                                            return 0;
                                        }
                                    } else {
                                        if (c.geti() == b.geti() && c.getj() == b.getj() && c.getk() == 0) {//(sppi|ppi?)
                                            switch (d.getL()) {
                                                case 0:
                                                    return sppisppideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                case 1:
                                                    if (d.getk() == 1) {
                                                        return sppippipzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                    } else {
                                                        return 0;
                                                    }
                                                default:
                                                    return 0;
                                            }
                                        } else {
                                            return 0;
                                        }
                                    }
                                default:
                                    System.err.println("oh no");
                                    return 0;
                            }
                        }
                    default:
                        System.err.println("oh no");
                        return 0;
                }

            case 1://(p?|??)
                if (a.getk() == 1) {//(pz?|??)
                    switch (b.getL()) {
                        case 0:
                            switch (c.getL()) {

                                case 0://(pzs|s?)

                                    switch (d.getL()) {

                                        case 0://(pzs|ss)
                                            return spzssderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);

                                        case 1:
                                            if (d.getk() == 1) {//(pzs|spz)
                                                return spzspzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                            } else {
                                                return 0;
                                            }
                                    }

                                case 1:
                                    if (c.getk() == 1) {//(pzs|pz?)

                                        switch (d.getL()) {

                                            case 0://(pzs|pzs)
                                                return spzspzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);

                                            case 1:
                                                if (d.getk() == 1) {//(pzs|pzpz)
                                                    return spzpzpzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                } else {//(pzs|pzppi) = 0
                                                    return 0;
                                                }
                                        }
                                    } else {//(pzs|ppi?)
                                        if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
                                            return spzppippideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                        } else {
                                            return 0;
                                        }
                                    }
                                default:
                                    System.err.println("oh no");
                                    return 0;
                            }
                        case 1:

                            if (b.getk() == 1) {//(pzpz|??)

                                switch (c.getL()) {

                                    case 0://(pzpz|s?)

                                        switch (d.getL()) {

                                            case 0://(pzpz|ss)
                                                return pzpzssderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);

                                            case 1:
                                                if (d.getk() == 1) {//(pzpz|spz)
                                                    return pzpzspzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                } else {
                                                    return 0;
                                                }
                                        }

                                    case 1:
                                        if (c.getk() == 1) {//(pzpz|pz?)

                                            switch (d.getL()) {

                                                case 0://(pzpz|pzs)
                                                    return pzpzspzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);

                                                case 1:
                                                    if (d.getk() == 1) {//(pzpz|pzpz)
                                                        return pzpzpzpzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                    } else {//(pzpz|pzppi) = 0
                                                        return 0;
                                                    }
                                            }
                                        } else {//(pzpz|ppi?)
                                            if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
                                                return pzpzppippideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                            } else {
                                                return 0;
                                            }
                                        }
                                    default:
                                        System.err.println("oh no");
                                        return 0;
                                }
                            } else {//(pzppi|??)

                                switch (c.getL()) {
                                    case 0://(pzppi|s?)
                                        if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(pzppi|sppi)
                                            return ppipzsppideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                        } else {
                                            return 0;
                                        }
                                    case 1:
                                        if (c.getk() == 1) {
                                            if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(pzppi|pzppi)
                                                return ppipzppipzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                            } else {
                                                return 0;
                                            }
                                        } else {
                                            if (c.geti() == b.geti() && c.getj() == b.getj() && c.getk() == 0) {//(pzppi|ppi?)
                                                switch (d.getL()) {
                                                    case 0:
                                                        return ppipzsppideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                    case 1:
                                                        if (d.getk() == 1) {
                                                            return ppipzppipzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                        } else {
                                                            return 0;
                                                        }
                                                    default:
                                                        return 0;
                                                }
                                            } else {
                                                return 0;
                                            }
                                        }
                                    default:
                                        System.err.println("oh no");
                                        return 0;
                                }
                            }
                    }
                } else {//(ppi?|??);

                    switch (b.getL()) {
                        case 0://(ppis|??)

                            switch (c.getL()) {
                                case 0://(ppis|s?)
                                    if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppis|sppi)
                                        return sppisppideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                    } else {
                                        return 0;
                                    }
                                case 1:
                                    if (c.getk() == 1) {
                                        if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppis|pzppi)
                                            return sppippipzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                        } else {
                                            return 0;
                                        }
                                    } else {
                                        if (c.geti() == a.geti() && c.getj() == a.getj() && c.getk() == 0) {//(ppis|ppi?)
                                            switch (d.getL()) {
                                                case 0:
                                                    return sppisppideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                case 1:
                                                    if (d.getk() == 1) {
                                                        return sppippipzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                    } else {
                                                        return 0;
                                                    }
                                                default:
                                                    return 0;
                                            }
                                        } else {
                                            return 0;
                                        }
                                    }
                                default:
                                    System.err.println("oh no");
                                    return 0;
                            }
                        case 1:
                            if (b.getk() == 1) {//(ppipz|??)
                                switch (c.getL()) {
                                    case 0://(ppipz|s?)
                                        if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppipz|sppi)
                                            return ppipzsppideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                        } else {
                                            return 0;
                                        }
                                    case 1:
                                        if (c.getk() == 1) {
                                            if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppipz|pzppi)
                                                return ppipzppipzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                            } else {
                                                return 0;
                                            }
                                        } else {
                                            if (c.geti() == a.geti() && c.getj() == a.getj() && c.getk() == 0) {//(ppipz|ppi?)
                                                switch (d.getL()) {
                                                    case 0:
                                                        return ppipzsppideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                    case 1:
                                                        if (d.getk() == 1) {
                                                            return ppipzppipzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                        } else {
                                                            return 0;
                                                        }
                                                    default:
                                                        return 0;
                                                }
                                            } else {
                                                return 0;
                                            }
                                        }
                                    default:
                                        System.err.println("oh no");
                                        return 0;
                                }

                            } else {

                                switch (c.getL()) {
                                    case 0://(ppippi|s?)
                                        switch (d.getL()) {
                                            case 0://(ppippi|ss)
                                                if (a.geti() == b.geti() && a.getj() == b.getj()) {
                                                    return ppippissderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                } else {
                                                    return 0;
                                                }
                                            case 1:
                                                if (d.getk() == 1 && a.geti() == b.geti() && a.getj() == b.getj()) {
                                                    return ppippispzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                } else {
                                                    return 0;
                                                }
                                        }

                                    case 1:
                                        if (c.getk() == 1) {
                                            switch (d.getL()) {
                                                case 0://(ppippi|pzs)
                                                    if (a.geti() == b.geti() && a.getj() == b.getj()) {
                                                        return ppippispzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                    } else {
                                                        return 0;
                                                    }

                                                case 1:
                                                    if (d.getk() == 1 && a.geti() == b.geti() && a.getj() == b.getj()) {//(ppippi|pzpz)
                                                        return ppippipzpzderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                    } else {
                                                        return 0;
                                                    }
                                            }
                                        } else {
                                            if (a.geti() == b.geti() && a.getj() == b.getj()) {//(pxpx|??) or (pypy|??)

                                                if (c.getL() == d.getL() && c.geti() == d.geti() && c.getj() == d.getj() && c.getk() == 0) {
                                                    if (a.geti() == c.geti()) {
                                                        return ppippippippideriv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                    } else {
                                                        return pxpxpypyderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                    }
                                                } else {
                                                    return 0;
                                                }

                                            } else {//(pxpy|??) or (pypx|??)
                                                if (c.getL() == d.getL() && c.geti() != d.geti() && c.getj() != d.getj() && c.getk() == 0) {
                                                    return pxpypxpyderiv2(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau1, tau2);
                                                }
                                            }
                                        }

                                }

                            }
                    }

                }
        }

        return 0;
    }

    public static double getGderiv2finite(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau1, int tau2) {
        double orig = NDDODerivative.getGderiv(a, b, c, d, tau1);

        double[] newcoords = a.getCoords().clone();

        newcoords[tau2] += 1E-6;

        NDDO6G anew = new NDDO6G(a, newcoords);
        NDDO6G bnew = new NDDO6G(b, newcoords);

        double perturbed = NDDODerivative.getGderiv(anew, bnew, c, d, tau1);

        return (perturbed - orig) / 1E-6;
    }

    public static double[] secondDerivativeDecompositionfinite(double[] point1, double[] point2, NDDO6G a, int tau1, int tau2) {
        if (a.getL() == 0) {
            return new double[]{0};
        }

        double[] orig = NDDODerivative.derivativeDecomposition(point1, point2, a, tau1);

        point1 = point1.clone();

        point1[tau2] += 1E-9;

        double[] perturbed = NDDODerivative.derivativeDecomposition(point1, point2, a, tau1);

        return new double[]{(perturbed[0] - orig[0]) / 1E-9, (perturbed[1] - orig[1]) / 1E-9, (perturbed[2] - orig[2]) / 1E-9};
    }

    public static double[] secondDerivativeDecomposition(double[] point1, double[] point2, NDDO6G a, int tau1, int tau2) {
        if (a.getL() == 0) {
            return new double[]{0};
        }

        int A = Math.min(tau1, tau2);
        int B = Math.max(tau1, tau2);

        double x = point2[0] - point1[0];
        double y = point2[1] - point1[1];
        double z = point2[2] - point1[2];

        double R = GTO.R(point1, point2);
        double Rxy = Math.sqrt(x * x + y * y);

        double[] returnval = new double[3];


        if (a.geti() == 1 && a.getL() == 1) {

            switch (A) {
                case 0:

                    switch (B) {
                        case 0://partial wrt x and x
                            returnval[0] = 3 * x * z / (R * R * R * Rxy) + 3 * x * z / (R * Rxy * Rxy * Rxy)
                                    - 2 * x * x * x * z / (R * R * R * Rxy * Rxy * Rxy)
                                    - 3 * x * x * x * z / (R * Rxy * Rxy * Rxy * Rxy * Rxy)
                                    - 3 * x * x * x * z / (R * R * R * R * R * Rxy);

                            returnval[1] = y / (Rxy * Rxy * Rxy) - 3 * x * x * y / (Rxy * Rxy * Rxy * Rxy * Rxy);

                            returnval[2] = 3 * x * x * x / (R * R * R * R * R) - 3 * x / (R * R * R);


                            return returnval;

                        case 1://partial wrt x and y
                            returnval[0] = z * y / (R * R * R * Rxy) + z * y / (R * Rxy * Rxy * Rxy)
                                    - 2 * x * x * y * z / (R * R * R * Rxy * Rxy * Rxy)
                                    - 3 * x * x * y * z / (R * R * R * R * R * Rxy)
                                    - 3 * x * x * y * z / (Rxy * Rxy * Rxy * Rxy * Rxy * R);

                            returnval[1] = x / (Rxy * Rxy * Rxy) - 3 * x * y * y / (Rxy * Rxy * Rxy * Rxy * Rxy);

                            returnval[2] = 3 * x * x * y / (R * R * R * R * R) - y / (R * R * R);


                            return returnval;

                        case 2://partial wrt x and z
                            returnval[0] = z * z / (R * R * R * Rxy) - 1 / (R * Rxy)
                                    + x * x / (R * R * R * Rxy) + x * x / (Rxy * Rxy * Rxy * R)
                                    - 3 * x * x * z * z / (R * R * R * R * R * Rxy)
                                    - x * x * z * z / (R * R * R * Rxy * Rxy * Rxy);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * x * z / (R * R * R * R * R) - z / (R * R * R);


                            return returnval;

                    }
                    break;
                case 1:

                    switch (B) {

                        case 1://partial wrt y and y
                            returnval[0] = x * z / (R * R * R * Rxy) + x * z / (R * Rxy * Rxy * Rxy)
                                    - 2 * x * y * y * z / (R * R * R * Rxy * Rxy * Rxy)
                                    - 3 * x * y * y * z / (R * R * R * R * R * Rxy)
                                    - 3 * x * y * y * z / (Rxy * Rxy * Rxy * Rxy * Rxy * R);

                            returnval[1] = 3 * y / (Rxy * Rxy * Rxy) - 3 * y * y * y / (Rxy * Rxy * Rxy * Rxy * Rxy);

                            returnval[2] = 3 * x * y * y / (R * R * R * R * R) - x / (R * R * R);

                            return returnval;

                        case 2://partial wrt y and z
                            returnval[0] = x * y / (R * R * R * Rxy) + x * y / (R * Rxy * Rxy * Rxy)
                                    - 3 * x * y * z * z / (R * R * R * R * R * Rxy)
                                    - x * y * z * z / (R * R * R * Rxy * Rxy * Rxy);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * y * z / (R * R * R * R * R);


                            return returnval;

                    }
                case 2:
                    returnval[0] = 3 * x * z / (R * R * R * Rxy) - 3 * x * z * z * z / (R * R * R * R * R * Rxy);

                    returnval[1] = 0;

                    returnval[2] = 3 * x * z * z / (R * R * R * R * R) - x / (R * R * R);


                    return returnval;
                default:
            }
        } else if (a.getj() == 1 && a.getL() == 1) {

            switch (A) {
                case 0:

                    switch (B) {
                        case 0: //partial wrt x and x
                            returnval[0] = y * z / (R * R * R * Rxy) + y * z / (R * Rxy * Rxy * Rxy)
                                    - 2 * x * x * y * z / (R * R * R * Rxy * Rxy * Rxy)
                                    - 3 * x * x * y * z / (R * R * R * R * R * Rxy)
                                    - 3 * x * x * y * z / (Rxy * Rxy * Rxy * Rxy * Rxy * R);

                            returnval[1] = -3 * x / (Rxy * Rxy * Rxy) + 3 * x * x * x / (Rxy * Rxy * Rxy * Rxy * Rxy);

                            returnval[2] = 3 * x * x * y / (R * R * R * R * R) - y / (R * R * R);

                            return returnval;

                        case 1://partial wrt x and y
                            returnval[0] = z * x / (R * R * R * Rxy) + z * x / (R * Rxy * Rxy * Rxy)
                                    - 2 * x * y * y * z / (R * R * R * Rxy * Rxy * Rxy)
                                    - 3 * x * y * y * z / (R * R * R * R * R * Rxy)
                                    - 3 * x * y * y * z / (Rxy * Rxy * Rxy * Rxy * Rxy * R);

                            returnval[1] = -y / (Rxy * Rxy * Rxy) + 3 * x * x * y / (Rxy * Rxy * Rxy * Rxy * Rxy);

                            returnval[2] = 3 * x * y * y / (R * R * R * R * R) - x / (R * R * R);


                            return returnval;

                        case 2://partial wrt x and z

                            returnval[0] = x * y / (R * R * R * Rxy) + x * y / (R * Rxy * Rxy * Rxy)
                                    - 3 * x * y * z * z / (R * R * R * R * R * Rxy)
                                    - x * y * z * z / (R * R * R * Rxy * Rxy * Rxy);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * y * z / (R * R * R * R * R);


                            return returnval;

                    }
                    break;
                case 1:

                    switch (B) {

                        case 1: //partial wrt y and y
                            returnval[0] = 3 * y * z / (R * R * R * Rxy) + 3 * y * z / (R * Rxy * Rxy * Rxy)
                                    - 2 * y * y * y * z / (R * R * R * Rxy * Rxy * Rxy)
                                    - 3 * y * y * y * z / (R * Rxy * Rxy * Rxy * Rxy * Rxy)
                                    - 3 * y * y * y * z / (R * R * R * R * R * Rxy);

                            returnval[1] = -x / (Rxy * Rxy * Rxy) + 3 * x * y * y / (Rxy * Rxy * Rxy * Rxy * Rxy);

                            returnval[2] = 3 * y * y * y / (R * R * R * R * R) - 3 * y / (R * R * R);

                            return returnval;

                        case 2://partial wrt y and z
                            returnval[0] = z * z / (R * R * R * Rxy) - 1 / (R * Rxy)
                                    + y * y / (R * R * R * Rxy) + y * y / (Rxy * Rxy * Rxy * R)
                                    - 3 * y * y * z * z / (R * R * R * R * R * Rxy)
                                    - y * y * z * z / (R * R * R * Rxy * Rxy * Rxy);

                            returnval[1] = 0;

                            returnval[2] = 3 * y * y * z / (R * R * R * R * R) - z / (R * R * R);


                            return returnval;

                    }
                case 2:
                    returnval[0] = 3 * y * z / (R * R * R * Rxy) - 3 * y * z * z * z / (R * R * R * R * R * Rxy);

                    returnval[1] = 0;

                    returnval[2] = 3 * y * z * z / (R * R * R * R * R) - y / (R * R * R);


                    return returnval;
                default:
            }
        } else if (a.getk() == 1 && a.getL() == 1) {

            switch (A) {
                case 0:

                    switch (B) {

                        case 0: //partial wrt x and x
                            returnval[0] = 3 * Rxy * x * x / (R * R * R * R * R) - Rxy / (R * R * R) + 1 / (Rxy * R)
                                    - 2 * x * x / (Rxy * R * R * R) - x * x / (R * Rxy * Rxy * Rxy);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * x * z / (R * R * R * R * R) - z / (R * R * R);

                            return returnval;

                        case 1://partial wrt x and y
                            returnval[0] = 3 * Rxy * x * y / (R * R * R * R * R) - 2 * x * y / (Rxy * R * R * R)
                                    - x * y / (R * Rxy * Rxy * Rxy);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * y * z / (R * R * R * R * R);


                            return returnval;

                        case 2://partial wrt x and z

                            returnval[0] = 3 * Rxy * x * z / (R * R * R * R * R) - x * z / (Rxy * R * R * R);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * z * z / (R * R * R * R * R) - x / (R * R * R);


                            return returnval;

                    }
                    break;
                case 1:

                    switch (B) {

                        case 1: //partial wrt y and y
                            returnval[0] = 3 * Rxy * y * y / (R * R * R * R * R) - Rxy / (R * R * R) + 1 / (Rxy * R)
                                    - 2 * y * y / (Rxy * R * R * R) - y * y / (R * Rxy * Rxy * Rxy);

                            returnval[1] = 0;

                            returnval[2] = 3 * y * y * z / (R * R * R * R * R) - z / (R * R * R);

                            return returnval;

                        case 2://partial wrt y and z

                            returnval[0] = 3 * Rxy * y * z / (R * R * R * R * R) - y * z / (Rxy * R * R * R);

                            returnval[1] = 0;

                            returnval[2] = 3 * y * z * z / (R * R * R * R * R) - y / (R * R * R);


                            return returnval;

                    }
                case 2:
                    returnval[0] = 3 * z * z * Rxy / (R * R * R * R * R) - Rxy / (R * R * R);

                    returnval[1] = 0;

                    returnval[2] = 3 * z * z * z / (R * R * R * R * R) - 3 * z / (R * R * R);

                    return returnval;
                default:
            }
        }

        System.err.println("oh no!");

        return new double[]{0, 0, 0};


    }

    public static double[] secondDerivativeDecomposition2finite(double[] point1, double[] point2, NDDO6G a, int tau1, int tau2) {
        if (a.getL() == 0) {
            return new double[]{0};
        }

        double[] orig = NDDODerivative.derivativeDecomposition2(point1, point2, a, tau1);

        point1 = point1.clone();

        point1[tau2] += 1E-9;

        double[] perturbed = NDDODerivative.derivativeDecomposition2(point1, point2, a, tau1);

        return new double[]{(perturbed[0] - orig[0]) / 1E-9, (perturbed[1] - orig[1]) / 1E-9, (perturbed[2] - orig[2]) / 1E-9};
    }

    public static double[] secondDerivativeDecomposition2(double[] point1, double[] point2, NDDO6G a, int tau1, int tau2) {
        if (a.getL() == 0) {
            return new double[]{0};
        }

        int A = Math.min(tau1, tau2);
        int B = Math.max(tau1, tau2);

        double x = point2[0] - point1[0];
        double y = point2[1] - point1[1];
        double z = point2[2] - point1[2];

        double R = GTO.R(point1, point2);
        double Rxz = Math.sqrt(x * x + z * z);

        double[] returnval = new double[3];

        if (a.geti() == 1 && a.getL() == 1) {

            switch (A) {
                case 0:

                    switch (B) {
                        case 0://partial wrt x and x
                            returnval[0] = 3 * x * x * x * y / (R * Rxz * Rxz * Rxz * Rxz * Rxz)
                                    + 3 * x * x * x * y / (R * R * R * R * R * Rxz)
                                    + 2 * x * x * x * y / (R * R * R * Rxz * Rxz * Rxz)
                                    - 3 * x * y / (R * R * R * Rxz)
                                    - 3 * x * y / (R * Rxz * Rxz * Rxz);

                            returnval[1] = z / (Rxz * Rxz * Rxz) - 3 * x * x * z / (Rxz * Rxz * Rxz * Rxz * Rxz);

                            returnval[2] = 3 * x * x * x / (R * R * R * R * R) - 3 * x / (R * R * R);


                            return returnval;

                        case 1://partial wrt x and y
                            returnval[0] = 1 / (R * Rxz) + 3 * x * x * y * y / (R * R * R * R * R * Rxz)
                                    + x * x * y * y / (R * R * R * Rxz * Rxz * Rxz)
                                    - x * x / (R * Rxz * Rxz * Rxz) - x * x / (R * R * R * Rxz)
                                    - y * y / (R * R * R * Rxz);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * x * y / (R * R * R * R * R) - y / (R * R * R);


                            return returnval;

                        case 2://partial wrt x and z
                            returnval[0] = 3 * x * x * y * z / (R * R * R * R * R * Rxz)
                                    + 2 * x * x * y * z / (R * R * R * Rxz * Rxz * Rxz)
                                    + 3 * x * x * y * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz)
                                    - y * z / (R * R * R * Rxz)
                                    - y * z / (R * Rxz * Rxz * Rxz);


                            returnval[1] = x / (Rxz * Rxz * Rxz) - 3 * x * z * z / (Rxz * Rxz * Rxz * Rxz * Rxz);

                            returnval[2] = 3 * x * x * z / (R * R * R * R * R) - z / (R * R * R);


                            return returnval;

                    }
                    break;
                case 1:

                    switch (B) {

                        case 1://partial wrt y and y
                            returnval[0] = 3 * x * y * y * y / (R * R * R * R * R * Rxz)
                                    - 3 * x * y / (R * R * R * Rxz);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * y * y / (R * R * R * R * R) - x / (R * R * R);

                            return returnval;

                        case 2://partial wrt y and z
                            returnval[0] = 3 * x * y * y * z / (R * R * R * R * R * Rxz)
                                    + x * y * y * z / (R * R * R * Rxz * Rxz * Rxz)
                                    - x * z / (R * Rxz * Rxz * Rxz)
                                    - x * z / (R * R * R * Rxz);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * y * z / (R * R * R * R * R);


                            return returnval;

                    }
                case 2:
                    returnval[0] = 3 * x * y * z * z / (R * R * R * R * R * Rxz)
                            + 3 * x * y * z * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz)
                            + 2 * x * y * z * z / (R * R * R * Rxz * Rxz * Rxz)
                            - x * y / (R * Rxz * Rxz * Rxz) - x * y / (R * R * R * Rxz);

                    returnval[1] = 3 * z / (Rxz * Rxz * Rxz) - 3 * z * z * z / (Rxz * Rxz * Rxz * Rxz * Rxz);

                    returnval[2] = 3 * x * z * z / (R * R * R * R * R) - x / (R * R * R);


                    return returnval;
                default:
            }
        } else if (a.getj() == 1 && a.getL() == 1) {

            switch (A) {
                case 0:

                    switch (B) {
                        case 0: //partial wrt x and x
                            returnval[0] = 2 * x * x / (R * R * R * Rxz) + x * x / (R * Rxz * Rxz * Rxz)
                                    + Rxz / (R * R * R) - 1 / (R * Rxz) - 3 * x * x * Rxz / (R * R * R * R * R);
                            returnval[1] = 0;

                            returnval[2] = 3 * x * x * y / (R * R * R * R * R) - y / (R * R * R);

                            return returnval;

                        case 1://partial wrt x and y
                            returnval[0] = x * y / (R * R * R * Rxz) - 3 * x * y * Rxz / (R * R * R * R * R);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * y * y / (R * R * R * R * R) - x / (R * R * R);


                            return returnval;

                        case 2://partial wrt x and z

                            returnval[0] = 2 * x * z / (R * R * R * Rxz) + x * z / (R * Rxz * Rxz * Rxz)
                                    - 3 * x * z * Rxz / (R * R * R * R * R);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * y * z / (R * R * R * R * R);


                            return returnval;

                    }
                    break;
                case 1:

                    switch (B) {

                        case 1: //partial wrt y and y
                            returnval[0] = Rxz / (R * R * R) - 3 * y * y * Rxz / (R * R * R * R * R);

                            returnval[1] = 0;

                            returnval[2] = 3 * y * y * y / (R * R * R * R * R) - 3 * y / (R * R * R);

                            return returnval;

                        case 2://partial wrt y and z
                            returnval[0] = y * z / (R * R * R * Rxz) - 3 * y * z * Rxz / (R * R * R * R * R);

                            returnval[1] = 0;

                            returnval[2] = 3 * y * y * z / (R * R * R * R * R) - z / (R * R * R);


                            return returnval;

                    }
                case 2:
                    returnval[0] = 2 * z * z / (R * R * R * Rxz) + z * z / (R * Rxz * Rxz * Rxz)
                            + Rxz / (R * R * R) - 1 / (R * Rxz) - 3 * z * z * Rxz / (R * R * R * R * R);

                    returnval[1] = 0;

                    returnval[2] = 3 * y * z * z / (R * R * R * R * R) - y / (R * R * R);


                    return returnval;
                default:
            }
        } else if (a.getk() == 1 && a.getL() == 1) {

            switch (A) {
                case 0:

                    switch (B) {

                        case 0: //partial wrt x and x
                            returnval[0] = 2 * x * x * y * z / (R * R * R * Rxz * Rxz * Rxz) + 3 * x * x * y * z / (R * R * R * R * R * Rxz)
                                    + 3 * x * x * y * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz) - y * z / (R * R * R * Rxz) - y * z / (R * Rxz * Rxz * Rxz);

                            returnval[1] = 3 * x * x * x / (Rxz * Rxz * Rxz * Rxz * Rxz) - 3 * x / (Rxz * Rxz * Rxz);

                            returnval[2] = 3 * x * x * z / (R * R * R * R * R) - z / (R * R * R);

                            return returnval;

                        case 1://partial wrt x and y
                            returnval[0] = 3 * x * y * y * z / (R * R * R * R * R * Rxz)
                                    + x * y * y * z / (R * R * R * Rxz * Rxz * Rxz)
                                    - x * z / (R * Rxz * Rxz * Rxz) - x * z / (R * R * R * Rxz);

                            returnval[1] = 0;

                            returnval[2] = 3 * x * y * z / (R * R * R * R * R);


                            return returnval;

                        case 2://partial wrt x and z

                            returnval[0] = 2 * x * y * z * z / (R * R * R * Rxz * Rxz * Rxz) + 3 * x * y * z * z / (R * R * R * R * R * Rxz)
                                    + 3 * x * y * z * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz) - x * y / (R * Rxz * Rxz * Rxz) - x * y / (R * R * R * Rxz);

                            returnval[1] = 3 * x * x * z / (Rxz * Rxz * Rxz * Rxz * Rxz) - z / (Rxz * Rxz * Rxz);

                            returnval[2] = 3 * x * z * z / (R * R * R * R * R) - x / (R * R * R);


                            return returnval;

                    }
                    break;
                case 1:

                    switch (B) {

                        case 1: //partial wrt y and y
                            returnval[0] = 3 * y * y * y * z / (R * R * R * R * R * Rxz) - 3 * y * z / (R * R * R * Rxz);

                            returnval[1] = 0;

                            returnval[2] = 3 * y * y * z / (R * R * R * R * R) - z / (R * R * R);

                            return returnval;

                        case 2://partial wrt y and z

                            returnval[0] = 1 / (R * Rxz) + 3 * y * y * z * z / (R * R * R * R * R * Rxz)
                                    + y * y * z * z / (R * R * R * Rxz * Rxz * Rxz) - y * y / (R * R * R * Rxz)
                                    - z * z / (R * Rxz * Rxz * Rxz) - z * z / (R * R * R * Rxz);

                            returnval[1] = 0;

                            returnval[2] = 3 * y * z * z / (R * R * R * R * R) - y / (R * R * R);


                            return returnval;

                    }
                case 2:
                    returnval[0] = 3 * y * z * z * z / (R * R * R * R * R * Rxz) + 2 * y * z * z * z / (R * R * R * Rxz * Rxz * Rxz)
                            + 3 * y * z * z * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz) - 3 * y * z / (R * Rxz * Rxz * Rxz) - 3 * y * z / (R * R * R * Rxz);

                    returnval[1] = 3 * x * z * z / (Rxz * Rxz * Rxz * Rxz * Rxz) - x / (Rxz * Rxz * Rxz);

                    returnval[2] = 3 * z * z * z / (R * R * R * R * R) - 3 * z / (R * R * R);

                    return returnval;
                default:
            }
        }

        System.err.println("oh no!");

        return new double[]{0, 0, 0};


    }

    public static double getGderiv2(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau1, int tau2) {


        double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
        double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
        double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
        double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

        double[] coeffAderiv1 = NDDODerivative.derivativeDecomposition(a.getCoords(), c.getCoords(), a, tau1);
        double[] coeffBderiv1 = NDDODerivative.derivativeDecomposition(a.getCoords(), c.getCoords(), b, tau1);
        double[] coeffCderiv1 = NDDODerivative.derivativeDecomposition(a.getCoords(), c.getCoords(), c, tau1);
        double[] coeffDderiv1 = NDDODerivative.derivativeDecomposition(a.getCoords(), c.getCoords(), d, tau1);

        double[] coeffAderiv2 = NDDODerivative.derivativeDecomposition(a.getCoords(), c.getCoords(), a, tau2);
        double[] coeffBderiv2 = NDDODerivative.derivativeDecomposition(a.getCoords(), c.getCoords(), b, tau2);
        double[] coeffCderiv2 = NDDODerivative.derivativeDecomposition(a.getCoords(), c.getCoords(), c, tau2);
        double[] coeffDderiv2 = NDDODerivative.derivativeDecomposition(a.getCoords(), c.getCoords(), d, tau2);

        double[] coeffAderiv = secondDerivativeDecomposition(a.getCoords(), c.getCoords(), a, tau1, tau2);
        double[] coeffBderiv = secondDerivativeDecomposition(a.getCoords(), c.getCoords(), b, tau1, tau2);
        double[] coeffCderiv = secondDerivativeDecomposition(a.getCoords(), c.getCoords(), c, tau1, tau2);
        double[] coeffDderiv = secondDerivativeDecomposition(a.getCoords(), c.getCoords(), d, tau1, tau2);

        if (Math.abs(a.getCoords()[0] - c.getCoords()[0]) < 1E-3 && Math.abs(a.getCoords()[1] - c.getCoords()[1]) < 1E-3) {

            coeffAderiv = secondDerivativeDecomposition2(a.getCoords(), c.getCoords(), a, tau1, tau2);
            coeffBderiv = secondDerivativeDecomposition2(a.getCoords(), c.getCoords(), b, tau1, tau2);
            coeffCderiv = secondDerivativeDecomposition2(a.getCoords(), c.getCoords(), c, tau1, tau2);
            coeffDderiv = secondDerivativeDecomposition2(a.getCoords(), c.getCoords(), d, tau1, tau2);

            coeffAderiv2 = NDDODerivative.derivativeDecomposition2(a.getCoords(), c.getCoords(), a, tau2);
            coeffBderiv2 = NDDODerivative.derivativeDecomposition2(a.getCoords(), c.getCoords(), b, tau2);
            coeffCderiv2 = NDDODerivative.derivativeDecomposition2(a.getCoords(), c.getCoords(), c, tau2);
            coeffDderiv2 = NDDODerivative.derivativeDecomposition2(a.getCoords(), c.getCoords(), d, tau2);

            coeffAderiv1 = NDDODerivative.derivativeDecomposition2(a.getCoords(), c.getCoords(), a, tau1);
            coeffBderiv1 = NDDODerivative.derivativeDecomposition2(a.getCoords(), c.getCoords(), b, tau1);
            coeffCderiv1 = NDDODerivative.derivativeDecomposition2(a.getCoords(), c.getCoords(), c, tau1);
            coeffDderiv1 = NDDODerivative.derivativeDecomposition2(a.getCoords(), c.getCoords(), d, tau1);

            coeffA = a.decomposition2(a.getCoords(), c.getCoords());
            coeffB = b.decomposition2(a.getCoords(), c.getCoords());
            coeffC = c.decomposition2(a.getCoords(), c.getCoords());
            coeffD = d.decomposition2(a.getCoords(), c.getCoords());


        }


        NDDO6G[] A = a.orbitalArray();
        NDDO6G[] B = b.orbitalArray();
        NDDO6G[] C = c.orbitalArray();
        NDDO6G[] D = d.orbitalArray();

        double sum = 0;

        for (int i = 0; i < coeffA.length; i++) {
            for (int j = 0; j < coeffB.length; j++) {
                for (int k = 0; k < coeffC.length; k++) {
                    for (int l = 0; l < coeffD.length; l++) {

                        double eri = NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k], D[l]);

                        double erideriv1 = NDDODerivative.LocalTwoCenterERIderiv(A[i], B[j], C[k], D[l], tau1);

                        double erideriv2 = NDDODerivative.LocalTwoCenterERIderiv(A[i], B[j], C[k], D[l], tau2);

                        double erideriv = LocalTwoCenterERIderiv2(A[i], B[j], C[k], D[l], tau1, tau2);

                        sum += coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] * eri * 27.21;
                        sum += coeffAderiv1[i] * coeffBderiv2[j] * coeffC[k] * coeffD[l] * eri * 27.21;
                        sum += coeffAderiv1[i] * coeffB[j] * coeffCderiv2[k] * coeffD[l] * eri * 27.21;
                        sum += coeffAderiv1[i] * coeffB[j] * coeffC[k] * coeffDderiv2[l] * eri * 27.21;
                        sum += coeffAderiv1[i] * coeffB[j] * coeffC[k] * coeffD[l] * erideriv2 * 27.21;

                        sum += coeffAderiv2[i] * coeffBderiv1[j] * coeffC[k] * coeffD[l] * eri * 27.21;
                        sum += coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] * eri * 27.21;
                        sum += coeffA[i] * coeffBderiv1[j] * coeffCderiv2[k] * coeffD[l] * eri * 27.21;
                        sum += coeffA[i] * coeffBderiv1[j] * coeffC[k] * coeffDderiv2[l] * eri * 27.21;
                        sum += coeffA[i] * coeffBderiv1[j] * coeffC[k] * coeffD[l] * erideriv2 * 27.21;

                        sum += coeffAderiv2[i] * coeffB[j] * coeffCderiv1[k] * coeffD[l] * eri * 27.21;
                        sum += coeffA[i] * coeffBderiv2[j] * coeffCderiv1[k] * coeffD[l] * eri * 27.21;
                        sum += coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] * eri * 27.21;
                        sum += coeffA[i] * coeffB[j] * coeffCderiv1[k] * coeffDderiv2[l] * eri * 27.21;
                        sum += coeffA[i] * coeffB[j] * coeffCderiv1[k] * coeffD[l] * erideriv2 * 27.21;

                        sum += coeffAderiv2[i] * coeffB[j] * coeffC[k] * coeffDderiv1[l] * eri * 27.21;
                        sum += coeffA[i] * coeffBderiv2[j] * coeffC[k] * coeffDderiv1[l] * eri * 27.21;
                        sum += coeffA[i] * coeffB[j] * coeffCderiv2[k] * coeffDderiv1[l] * eri * 27.21;
                        sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] * eri * 27.21;
                        sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv1[l] * erideriv2 * 27.21;

                        sum += coeffAderiv2[i] * coeffB[j] * coeffC[k] * coeffD[l] * erideriv1 * 27.21;
                        sum += coeffA[i] * coeffBderiv2[j] * coeffC[k] * coeffD[l] * erideriv1 * 27.21;
                        sum += coeffA[i] * coeffB[j] * coeffCderiv2[k] * coeffD[l] * erideriv1 * 27.21;
                        sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv2[l] * erideriv1 * 27.21;
                        sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] * erideriv * 27.21;


                    }
                }
            }
        }

        return sum;
    }

    private static double Ederiv2(int atomnum1, int atomnum2, int[][] index, DoubleMatrix densityMatrix, NDDOAtom[] atoms, NDDO6G[] orbitals, int tau1, int tau2) {

        double e = 0;

        for (int i : index[atomnum1]) {
            for (int j : index[atomnum1]) {
                if (i != -1 && j != -1) {
                    e += densityMatrix.get(i, j) * atoms[atomnum2].Vderiv2(orbitals[i], orbitals[j], tau1, tau2);
                }
            }
        }

        for (int k : index[atomnum2]) {
            for (int l : index[atomnum2]) {
                if (k != -1 && l != -1) {
                    e += densityMatrix.get(k, l) * atoms[atomnum1].Vderiv2(orbitals[k], orbitals[l], tau1, tau2);
                }
            }
        }

        for (int i : index[atomnum1]) {
            for (int k : index[atomnum2]) {
                if (i != -1 && k != -1) {
                    e += 2 * densityMatrix.get(i, k) * NDDO6G.betaderiv2(orbitals[i], orbitals[k], tau1, tau2);
                }
            }
        }

        for (int i : index[atomnum1]) {
            for (int j : index[atomnum1]) {
                for (int k : index[atomnum2]) {
                    for (int l : index[atomnum2]) {
                        if (i != -1 && j != -1 && k != -1 && l != -1) {
                            e += (densityMatrix.get(i, j) * densityMatrix.get(k, l) - densityMatrix.get(i, k) * 0.5 * densityMatrix.get(j, l))
                                    * NDDOSecondDerivative.getGderiv2(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau1, tau2);
                        }
                    }
                }
            }
        }

        return e;


    }

    private static double Ederiv2(int atomnum1, int atomnum2, int[][] index, DoubleMatrix alphaDensity, DoubleMatrix betaDensity, NDDOAtom[] atoms, NDDO6G[] orbitals, int tau1, int tau2) {

        double e = 0;

        for (int i : index[atomnum1]) {
            for (int j : index[atomnum1]) {
                if (i != -1 && j != -1) {
                    e += (alphaDensity.get(i, j) + betaDensity.get(i, j)) * atoms[atomnum2].Vderiv2(orbitals[i], orbitals[j], tau1, tau2);
                }
            }
        }

        for (int k : index[atomnum2]) {
            for (int l : index[atomnum2]) {
                if (k != -1 && l != -1) {
                    e += (alphaDensity.get(k, l) + betaDensity.get(k, l)) * atoms[atomnum1].Vderiv2(orbitals[k], orbitals[l], tau1, tau2);
                }
            }
        }

        for (int i : index[atomnum1]) {
            for (int k : index[atomnum2]) {
                if (i != -1 && k != -1) {
                    e += 2 * (alphaDensity.get(i, k) + betaDensity.get(i, k)) * NDDO6G.betaderiv2(orbitals[i], orbitals[k], tau1, tau2);
                }
            }
        }

        for (int i : index[atomnum1]) {
            for (int j : index[atomnum1]) {
                for (int k : index[atomnum2]) {
                    for (int l : index[atomnum2]) {
                        if (i != -1 && j != -1 && k != -1 && l != -1) {
                            e += ((alphaDensity.get(i, j) + betaDensity.get(i, j)) * (alphaDensity.get(k, l) + betaDensity.get(k, l)) - alphaDensity.get(i, k) * alphaDensity.get(j, l) - betaDensity.get(i, k) * betaDensity.get(j, l))
                                    * NDDOSecondDerivative.getGderiv2(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau1, tau2);
                        }
                    }
                }
            }
        }

        return e;


    }


    public static double hessianfinite(NDDOAtom[] atoms, NDDOSolutionRestricted soln, int atomnum1, int tau1, int atomnum2, int tau2) {

        double initval = NDDODerivative.gradient(atoms, soln, atomnum1, tau1);

        NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum2, tau2);

        NDDOSolutionRestricted newsoln = new NDDOSolutionRestricted(newatoms, soln.charge);

        double finalval = NDDODerivative.gradient(newatoms, newsoln, atomnum1, tau1);

        return 1E7 * (finalval - initval);


    }

    public static double hessianfinite(NDDOAtom[] atoms, NDDOSolutionUnrestricted soln, int atomnum1, int tau1, int atomnum2, int tau2) {

        double initval = NDDODerivative.gradientUnrestricted(atoms, soln, atomnum1, tau1);

        NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum2, tau2);

        NDDOSolutionUnrestricted newsoln = new NDDOSolutionUnrestricted(newatoms, soln.charge, soln.multiplicity);

        double finalval = NDDODerivative.gradientUnrestricted(newatoms, newsoln, atomnum1, tau1);

        return 1E7 * (finalval - initval);


    }

    public static DoubleMatrix hessianroutine(NDDOAtom[] atoms, NDDOSolutionRestricted soln, DoubleMatrix[] fockderivstatic) {

        DoubleMatrix[] densityderivs = new DoubleMatrix[fockderivstatic.length];

//        DoubleMatrix A = DoubleMatrix.zeros(NOcc * NVirt, NOcc * NVirt);
//
//        DoubleMatrix F = DoubleMatrix.zeros(NOcc * NVirt, fockderivstatic.length);
//
//        int count1 = 0;
//
//        for (int i = 0; i < NOcc; i++) {
//            for (int j = 0; j < NVirt; j++) {
//
//                int count2 = 0;
//
//                double energydiff = soln.E.get(NOcc + j) - soln.E.get(i);
//
//
//
//                for (int index = 0; index < fockderivstatic.length; index++) {
//
//                    double element = 0;
//
//                    for (int u = 0; u < soln.orbitals.length; u++) {
//                        for (int v = 0; v < soln.orbitals.length; v++) {
//                            element += soln.C.get(i, u) * soln.C.get(j + NOcc, v) * fockderivstatic[index].get(u, v);
//                        }
//                    }
//
//
//                    F.put(count1, index, element);
//                }
//
//
//                for (int k = 0; k < NOcc; k++) {
//                    for (int l = 0; l < NVirt; l++) {
//
//                        double coeff = 4 * ERIMOBasis(soln, i, j + NOcc, k, l + NOcc) - ERIMOBasis(soln, i, k, j + NOcc, l + NOcc) - ERIMOBasis(soln, i, l + NOcc, j + NOcc, k);
//
//                        if (i == k && j == l) {
//                            coeff += energydiff;
//                        }
//
//                        A.put(count1, count2, coeff);
//
//                        count2++;
//                    }
//                    //System.out.println ("A row of A is done");
//                }
//
//
//
//
//                count1++;
//            }
//        }
//
//        //System.out.println ("A is done");
//
//
//        if (NOcc == 0) {
//            for (int index = 0; index < densityderivs.length; index++) {
//                densityderivs[index] = DoubleMatrix.zeros (soln.orbitals.length, soln.orbitals.length);
//            }
//        }
//        else {
//            DoubleMatrix sol = Solve.solve(A, F);
//
//            System.err.println (sol);
//
//            //System.err.println ("This is how slow it is, you idiot");
//
//
//
//            for (int index = 0; index < densityderivs.length; index++) {
//                DoubleMatrix  densityderiv = DoubleMatrix.zeros (soln.orbitals.length, soln.orbitals.length);
//
//                for (int u = 0; u < densityderiv.rows; u++) {
//                    for (int v = 0; v < densityderiv.columns; v++) {
//                        double sum = 0;
//                        int count = 0;
//                        for (int i = 0; i < NOcc; i++) {
//                            for (int j = 0; j < NVirt; j++) {
//                                sum -= 2 * (soln.C.get(i, u) * soln.C.get(j + NOcc, v) + soln.C.get(j + NOcc, u) * soln.C.get(i, v)) * sol.get(count, index);
//                                count++;
//                            }
//                        }
//
//                        densityderiv.put(u, v, sum);
//                    }
//                }
//
//                densityderivs[index] = densityderiv;
//            }
//        }


        int count = 0;

        for (int a = 0; a < atoms.length; a++) {
            for (int tau = 0; tau < 3; tau++) {
                densityderivs[count] = NDDODerivative.densitymatrixderivfinite(atoms, soln, a, tau);
                densityderiv(atoms, soln, a, fockderivstatic[3 * a + tau]);
                count++;
            }
        }

        DoubleMatrix hessian = new DoubleMatrix(densityderivs.length, densityderivs.length);

        for (int i = 0; i < hessian.rows; i++) {
            for (int j = i; j < hessian.rows; j++) {

                double E = 0;

                int atomnum1 = i / 3;

                int atomnum2 = j / 3;

                int tau1 = i - 3 * atomnum1;

                int tau2 = j - 3 * atomnum2;

                if (atomnum1 == atomnum2) {
                    for (int a = 0; a < atoms.length; a++) {
                        if (a != atomnum1) {
                            E += Ederiv2(atomnum1, a, soln.index, soln.densityMatrix(), atoms, soln.orbitals, tau1, tau2);
                            E += atoms[atomnum1].crfDeriv2(atoms[a], tau1, tau2);
                        }
                    }
                } else {
                    E = -Ederiv2(atomnum1, atomnum2, soln.index, soln.densityMatrix(), atoms, soln.orbitals, tau1, tau2) - atoms[atomnum1].crfDeriv2(atoms[atomnum2], tau1, tau2);

                }


                for (int I = 0; I < soln.orbitals.length; I++) {
                    for (int J = 0; J < soln.orbitals.length; J++) {
                        E += fockderivstatic[i].get(I, J) * densityderivs[j].get(I, J);
                    }
                }

                hessian.put(i, j, E);
                hessian.put(j, i, E);
            }
        }

        return hessian;


    }

    public static DoubleMatrix densityderiv(NDDOAtom[] atoms, NDDOSolutionRestricted soln, int a, DoubleMatrix fockderivstatic) {

        int NOcc = (int) (soln.nElectrons / 2.0);

        int NVirt = soln.orbitals.length - NOcc;

        DoubleMatrix F = DoubleMatrix.zeros(NOcc * NVirt, 1);

        int count1 = 0;

        for (int i = 0; i < NOcc; i++) {
            for (int j = 0; j < NVirt; j++) {

                double element = 0;

                for (int u = 0; u < soln.orbitals.length; u++) {
                    for (int v = 0; v < soln.orbitals.length; v++) {
                        element += soln.C.get(i, u) * soln.C.get(j + NOcc, v) * fockderivstatic.get(u, v);
                    }
                }


                F.put(count1, 0, element);

                count1++;
            }
        }

        DoubleMatrix d = F.dup();

        if (d.rows == 0) {
            return new DoubleMatrix(0, 0);
        }

        DoubleMatrix mat = Singular.fullSVD(d)[0];

        ArrayList<DoubleMatrix> searchdirs = new ArrayList<>();

        for (int i = 0; i < mat.columns; i++) {

            DoubleMatrix vec = mat.getColumn(i);

            if (mag(vec) > 0) {
                searchdirs.add(vec);
            }
        }

        ArrayList<DoubleMatrix> responses = new ArrayList<>();

        for (int i = 0; i < searchdirs.size(); i++) {
            responses.add (computeResponseVectors(searchdirs.get(i), soln));
        }



        return null;


    }

    private static DoubleMatrix computeResponseVectors(DoubleMatrix x, NDDOSolutionRestricted soln) {

        int NOcc = (int) (soln.nElectrons / 2.0);

        int NVirt = soln.orbitals.length - NOcc;

        DoubleMatrix densityMatrixDeriv = DoubleMatrix.zeros(soln.orbitals.length, soln.orbitals.length);

        for (int u = 0; u < densityMatrixDeriv.rows; u++) {
            for (int v = 0; v < densityMatrixDeriv.columns; v++) {
                double sum = 0;
                int count = 0;
                for (int i = 0; i < NOcc; i++) {
                    for (int j = 0; j < NVirt; j++) {
                        sum -= 2 * (soln.C.get(i, u) * soln.C.get(j + NOcc, v) + soln.C.get(j + NOcc, u) * soln.C.get(i, v)) * x.get(count, 0);
                        count++;
                    }
                }

                densityMatrixDeriv.put(u, v, sum);
            }
        }

        DoubleMatrix responsematrix = DoubleMatrix.zeros(soln.orbitals.length, soln.orbitals.length);

        double[] integralArray = soln.integralArray;

        int integralcount = 0;

        for (int j = 0; j < soln.orbitals.length; j++) {
            for (int k = j; k < soln.orbitals.length; k++) {
                double val = 0;
                if (j == k) {

                    for (int l : soln.index[soln.atomNumber[j]]) {
                        if (l > -1) {
                            val += densityMatrixDeriv.get(l, l) * integralArray[integralcount];
                            integralcount++;
                        }
                    }

                    for (int l : soln.missingIndex[soln.atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : soln.missingIndex[soln.atomNumber[j]]) {
                                if (m > -1) {
                                    if (soln.atomNumber[l] == soln.atomNumber[m]) {
                                        val += densityMatrixDeriv.get(l, m) * integralArray[integralcount];
                                        integralcount++;
                                    }
                                }

                            }
                        }
                    }
                } else if (soln.atomNumber[j] == soln.atomNumber[k]) {
                    val += densityMatrixDeriv.get(j, k) * integralArray[integralcount];
                    integralcount++;

                    for (int l : soln.missingIndex[soln.atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : soln.missingIndex[soln.atomNumber[j]]) {
                                if (m > -1) {
                                    if (soln.atomNumber[l] == soln.atomNumber[m]) {
                                        val += densityMatrixDeriv.get(l, m) * integralArray[integralcount];
                                        integralcount++;
                                    }
                                }

                            }
                        }
                    }
                } else {
                    for (int l : soln.index[soln.atomNumber[j]]) {
                        if (l > -1) {
                            for (int m : soln.index[soln.atomNumber[k]]) {
                                if (m > -1) {
                                    val += densityMatrixDeriv.get(l, m) * integralArray[integralcount];
                                    integralcount++;
                                }
                            }
                        }
                    }
                }

                responsematrix.put(j, k, val);
                responsematrix.put(k, j, val);
            }
        }

        DoubleMatrix R = DoubleMatrix.zeros(NOcc * NVirt, 1);

        int count1 = 0;

        for (int i = 0; i < NOcc; i++) {
            for (int j = 0; j < NVirt; j++) {

                double element = 0;

                for (int u = 0; u < soln.orbitals.length; u++) {
                    for (int v = 0; v < soln.orbitals.length; v++) {
                        element += soln.C.get(i, u) * soln.C.get(j + NOcc, v) * responsematrix.get(u, v);
                    }
                }


                R.put(count1, 0, element);

                count1++;
            }
        }

        DoubleMatrix p = new DoubleMatrix (NOcc * NVirt, 1);

        int counter = 0;

        for (int i = 0; i < NOcc; i++) {
            for (int j = 0; j < NVirt; j++) {
                p.put (counter, 0, -R.get(counter, 0) + (soln.E.get(j + NOcc) - soln.E.get(i)) * x.get(counter));
                counter++;
            }
        }


        return p;
    }

    private static double mag(DoubleMatrix gradient) {

        double sum = 0;
        for (int i = 0; i < gradient.length; i++) {
            sum += gradient.get(i) * gradient.get(i);
        }

        return Math.sqrt(sum);
    }

    public static DoubleMatrix hessianroutine(NDDOAtom[] atoms, NDDOSolutionUnrestricted soln, DoubleMatrix[] fockderivstaticalpha, DoubleMatrix[] fockderivstaticbeta) {


        DoubleMatrix[] densityderivsalpha = new DoubleMatrix[fockderivstaticalpha.length];
        DoubleMatrix[] densityderivsbeta = new DoubleMatrix[fockderivstaticbeta.length];


        int count = 0;

        for (int a = 0; a < atoms.length; a++) {
            for (int tau = 0; tau < 3; tau++) {

                DoubleMatrix[] matrices = NDDODerivative.densitymatrixderivfinite(atoms, soln, a, tau);
                densityderivsalpha[count] = matrices[0];
                densityderivsbeta[count] = matrices[1];
                count++;
            }
        }

        DoubleMatrix hessian = new DoubleMatrix(densityderivsalpha.length, densityderivsalpha.length);

        for (int i = 0; i < hessian.rows; i++) {
            for (int j = i; j < hessian.rows; j++) {

                double E = 0;

                int atomnum1 = i / 3;

                int atomnum2 = j / 3;

                int tau1 = i - 3 * atomnum1;

                int tau2 = j - 3 * atomnum2;

                if (atomnum1 == atomnum2) {
                    for (int a = 0; a < atoms.length; a++) {
                        if (a != atomnum1) {
                            E += Ederiv2(atomnum1, a, soln.index, soln.alphaDensity(), soln.betaDensity(), atoms, soln.orbitals, tau1, tau2);
                            E += atoms[atomnum1].crfDeriv2(atoms[a], tau1, tau2);
                        }
                    }
                } else {
                    E = -Ederiv2(atomnum1, atomnum2, soln.index, soln.alphaDensity(), soln.betaDensity(), atoms, soln.orbitals, tau1, tau2) - atoms[atomnum1].crfDeriv2(atoms[atomnum2], tau1, tau2);

                }


                for (int I = 0; I < soln.orbitals.length; I++) {
                    for (int J = 0; J < soln.orbitals.length; J++) {
                        E += fockderivstaticalpha[i].get(I, J) * densityderivsalpha[j].get(I, J);
                        E += fockderivstaticbeta[i].get(I, J) * densityderivsbeta[j].get(I, J);
                    }
                }

                hessian.put(i, j, E);
                hessian.put(j, i, E);
            }
        }

        return hessian;


    }

    private static double ERIMOBasis(NDDOSolutionRestricted soln, int i, int j, int k, int l) {

        double eri = 0;

        for (int a = 0; a < soln.atoms.length; a++) {
            for (int b = 0; b < soln.atoms.length; b++) {
                for (int u : soln.index[a]) {
                    for (int v : soln.index[a]) {
                        for (int r : soln.index[b]) {
                            for (int s : soln.index[b]) {
                                if (u != -1 && v != -1 && r != -1 && s != -1) {
                                    if (a == b) {
                                        eri += soln.C.get(i, u) * soln.C.get(j, v) * soln.C.get(k, r) * soln.C.get(l, s) * NDDO6G.OneCenterERI(soln.orbitals[u], soln.orbitals[v], soln.orbitals[r], soln.orbitals[s]);
                                    } else {
                                        eri += soln.C.get(i, u) * soln.C.get(j, v) * soln.C.get(k, r) * soln.C.get(l, s) * NDDO6G.getG(soln.orbitals[u], soln.orbitals[v], soln.orbitals[r], soln.orbitals[s]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return eri;

    }


}
