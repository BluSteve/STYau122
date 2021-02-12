package nddoparam;

import org.jblas.DoubleMatrix;
import scf.GTO;
import scf.LCGTO;
import scf.Utils;

import java.util.Arrays;

public class NDDODerivative {

    private static double qqderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a00 = p01 + p02;
        return (xB[tau] - xA[tau]) * Math.pow(R * R + a00 * a00, -1.5);
    }

    private static double quzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a01 = p01 + p12;
        return 0.5 / R * (xB[tau] - xA[tau]) * (R + D12) * Math.pow((R + D12) * (R + D12) + a01 * a01, -1.5)
                - 0.5 / R * (xB[tau] - xA[tau]) * (R - D12) * Math.pow((R - D12) * (R - D12) + a01 * a01, -1.5);
    }

    private static double qQpipideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a02 = p01 + p22;
        return 0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D22 * D22 + a02 * a02, -1.5)
                - 0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + a02 * a02, -1.5);
    }

    private static double qQzzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a02 = p01 + p22;
        return 0.25 / R * (xB[tau] - xA[tau]) * (R + 2 * D22) * Math.pow((R + 2 * D22) * (R + 2 * D22) + a02 * a02, -1.5)
                + 0.25 / R * (xB[tau] - xA[tau]) * (R - 2 * D22) * Math.pow((R - 2 * D22) * (R - 2 * D22) + a02 * a02, -1.5)
                - 0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + a02 * a02, -1.5);
    }

    private static double upiupideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a11 = p11 + p12;
        return 0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + (D11 - D12) * (D11 - D12) + a11 * a11, -1.5)
                - 0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + (D11 + D12) * (D11 + D12) + a11 * a11, -1.5);
    }

    private static double uzuzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a11 = p11 + p12;
        return 0.25 / R * (xB[tau] - xA[tau]) * (R + D11 - D12) * Math.pow((R + D11 - D12) * (R + D11 - D12) + a11 * a11, -1.5)
                - 0.25 / R * (xB[tau] - xA[tau]) * (R + D11 + D12) * Math.pow((R + D11 + D12) * (R + D11 + D12) + a11 * a11, -1.5)
                - 0.25 / R * (xB[tau] - xA[tau]) * (R - D11 - D12) * Math.pow((R - D11 - D12) * (R - D11 - D12) + a11 * a11, -1.5)
                + 0.25 / R * (xB[tau] - xA[tau]) * (R - D11 + D12) * Math.pow((R - D11 + D12) * (R - D11 + D12) + a11 * a11, -1.5);
    }

    private static double upiQpizderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a12 = p11 + p22;
        return -0.25 / R * (xB[tau] - xA[tau]) * (R - D22) * Math.pow((R - D22) * (R - D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -1.5)
                + 0.25 / R * (xB[tau] - xA[tau]) * (R - D22) * Math.pow((R - D22) * (R - D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -1.5)
                + 0.25 / R * (xB[tau] - xA[tau]) * (R + D22) * Math.pow((R + D22) * (R + D22) + (D11 - D22) * (D11 - D22) + a12 * a12, -1.5)
                - 0.25 / R * (xB[tau] - xA[tau]) * (R + D22) * Math.pow((R + D22) * (R + D22) + (D11 + D22) * (D11 + D22) + a12 * a12, -1.5);
    }

    private static double uzQpipideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a12 = p11 + p22;
        return -0.25 / R * (xB[tau] - xA[tau]) * (R + D11) * Math.pow((R + D11) * (R + D11) + 4 * D22 * D22 + a12 * a12, -1.5)
                + 0.25 / R * (xB[tau] - xA[tau]) * (R - D11) * Math.pow((R - D11) * (R - D11) + 4 * D22 * D22 + a12 * a12, -1.5)
                + 0.25 / R * (xB[tau] - xA[tau]) * (R + D11) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
                - 0.25 / R * (xB[tau] - xA[tau]) * (R - D11) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
    }

    private static double uzQzzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a12 = p11 + p22;
        return -0.125 / R * (xB[tau] - xA[tau]) * (R + D11 - 2 * D22) * Math.pow((R + D11 - 2 * D22) * (R + D11 - 2 * D22) + a12 * a12, -1.5)
                + 0.125 / R * (xB[tau] - xA[tau]) * (R - D11 - 2 * D22) * Math.pow((R - D11 - 2 * D22) * (R - D11 - 2 * D22) + a12 * a12, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R + D11 + 2 * D22) * Math.pow((R + D11 + 2 * D22) * (R + D11 + 2 * D22) + a12 * a12, -1.5)
                + 0.125 / R * (xB[tau] - xA[tau]) * (R - D11 + 2 * D22) * Math.pow((R - D11 + 2 * D22) * (R - D11 + 2 * D22) + a12 * a12, -1.5)
                + 0.25 / R * (xB[tau] - xA[tau]) * (R + D11) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
                - 0.25 / R * (xB[tau] - xA[tau]) * (R - D11) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
    }

    private static double QpipiQpipideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a22 = p21 + p22;
        return 0.125 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
                + 0.125 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * (D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
                - 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
                - 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5)
                + 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
    }

    private static double QxxQyyderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a22 = p21 + p22;
        return 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, -1.5)
                - 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
                - 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5)
                + 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
    }

    private static double QpipiQzzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a22 = p21 + p22;
        return 0.125 / R * (xB[tau] - xA[tau]) * (R - 2 * D22) * Math.pow((R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5)
                + 0.125 / R * (xB[tau] - xA[tau]) * (R + 2 * D22) * Math.pow((R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R - 2 * D22) * Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R + 2 * D22) * Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5)
                - 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
                + 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
    }

    private static double QzzQzzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a22 = p21 + p22;
        return 0.0625 / R * (xB[tau] - xA[tau]) * (R + 2 * D21 - 2 * D22) * Math.pow((R + 2 * D21 - 2 * D22) * (R + 2 * D21 - 2 * D22) + a22 * a22, -1.5)
                + 0.0625 / R * (xB[tau] - xA[tau]) * (R + 2 * D21 + 2 * D22) * Math.pow((R + 2 * D21 + 2 * D22) * (R + 2 * D21 + 2 * D22) + a22 * a22, -1.5)
                + 0.0625 / R * (xB[tau] - xA[tau]) * (R - 2 * D21 - 2 * D22) * Math.pow((R - 2 * D21 - 2 * D22) * (R - 2 * D21 - 2 * D22) + a22 * a22, -1.5)
                + 0.0625 / R * (xB[tau] - xA[tau]) * (R - 2 * D21 + 2 * D22) * Math.pow((R - 2 * D21 + 2 * D22) * (R - 2 * D21 + 2 * D22) + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R + 2 * D21) * Math.pow((R + 2 * D21) * (R + 2 * D21) + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R - 2 * D21) * Math.pow((R - 2 * D21) * (R - 2 * D21) + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R + 2 * D22) * Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R - 2 * D22) * Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5)
                + 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
    }

    private static double QpizQpizderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {

        double R = GTO.R(xA, xB);
        double a22 = p21 + p22;
        return 0.125 / R * (xB[tau] - xA[tau]) * (R + D21 - D22) * Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R + D21 - D22) * Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R + D21 + D22) * Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
                + 0.125 / R * (xB[tau] - xA[tau]) * (R + D21 + D22) * Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R - D21 - D22) * Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
                + 0.125 / R * (xB[tau] - xA[tau]) * (R - D21 - D22) * Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
                + 0.125 / R * (xB[tau] - xA[tau]) * (R - D21 + D22) * Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
                - 0.125 / R * (xB[tau] - xA[tau]) * (R - D21 + D22) * Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5);
    }

    private static double ssssderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double ssppippideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double sspzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double ppippissderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
    }

    private static double pzpzssderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
    }

    private static double ppippippippideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QpipiQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double pxpxpypyderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QxxQyyderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double ppippipzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QpipiQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double pzpzppippideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QpipiQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
    }

    private static double pzpzpzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return qqderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QzzQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double spzssderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return -quzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
    }

    private static double spzppippideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return -quzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + uzQpipideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double spzpzpzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return -quzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + uzQzzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double ssspzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return quzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double ppippispzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return quzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) - uzQpipideriv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
    }

    private static double pzpzspzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return quzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) - uzQzzderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
    }

    private static double sppisppideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return upiupideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double spzspzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return uzuzderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double sppippipzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return upiQpizderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double ppipzsppideriv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return -upiQpizderiv(p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
    }

    private static double ppipzppipzderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return QpizQpizderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
    }

    private static double pxpypxpyderiv(double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
        return 0.5 * (ppippippippideriv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) - pxpxpypyderiv(p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau));
    }

    protected static double LocalTwoCenterERIderiv(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau) {

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
                                        return ssssderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);

                                    case 1:
                                        if (d.getk() == 1) {//(ss|spz)
                                            return ssspzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
                                            return ssspzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);

                                        case 1:
                                            if (d.getk() == 1) {//(ss|pzpz)
                                                return sspzpzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                            } else {//(ss|pzppi) = 0
                                                return 0;
                                            }
                                        default:
                                            return 0;
                                    }
                                } else {//(ss|ppi?)

                                    if (d.getL() == 1 && d.getk() == 0 && c.geti() == d.geti() && c.getj() == d.getj()) {//(ss|ppippi)
                                        return ssppippideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
                                            return spzssderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);

                                        case 1:
                                            if (d.getk() == 1) {//(spz|spz)
                                                return spzspzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                            } else {
                                                return 0;
                                            }
                                    }

                                case 1:
                                    if (c.getk() == 1) {//(spz|pz?)

                                        switch (d.getL()) {

                                            case 0://(spz|pzs)
                                                return spzspzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);

                                            case 1:
                                                if (d.getk() == 1) {//(spz|pzpz)
                                                    return spzpzpzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                } else {//(spz|pzppi) = 0
                                                    return 0;
                                                }
                                        }
                                    } else {//(spz|ppi?)
                                        if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
                                            return spzppippideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
                                        return sppisppideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                    } else {
                                        return 0;
                                    }
                                case 1:
                                    if (c.getk() == 1) {
                                        if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(sppi|pzppi)
                                            return sppippipzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                        } else {
                                            return 0;
                                        }
                                    } else {
                                        if (c.geti() == b.geti() && c.getj() == b.getj() && c.getk() == 0) {//(sppi|ppi?)
                                            switch (d.getL()) {
                                                case 0:
                                                    return sppisppideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                case 1:
                                                    if (d.getk() == 1) {
                                                        return sppippipzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
                                            return spzssderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);

                                        case 1:
                                            if (d.getk() == 1) {//(pzs|spz)
                                                return spzspzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                            } else {
                                                return 0;
                                            }
                                    }

                                case 1:
                                    if (c.getk() == 1) {//(pzs|pz?)

                                        switch (d.getL()) {

                                            case 0://(pzs|pzs)
                                                return spzspzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);

                                            case 1:
                                                if (d.getk() == 1) {//(pzs|pzpz)
                                                    return spzpzpzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                } else {//(pzs|pzppi) = 0
                                                    return 0;
                                                }
                                        }
                                    } else {//(pzs|ppi?)
                                        if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
                                            return spzppippideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
                                                return pzpzssderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);

                                            case 1:
                                                if (d.getk() == 1) {//(pzpz|spz)
                                                    return pzpzspzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                } else {
                                                    return 0;
                                                }
                                        }

                                    case 1:
                                        if (c.getk() == 1) {//(pzpz|pz?)

                                            switch (d.getL()) {

                                                case 0://(pzpz|pzs)
                                                    return pzpzspzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);

                                                case 1:
                                                    if (d.getk() == 1) {//(pzpz|pzpz)
                                                        return pzpzpzpzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                    } else {//(pzpz|pzppi) = 0
                                                        return 0;
                                                    }
                                            }
                                        } else {//(pzpz|ppi?)
                                            if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
                                                return pzpzppippideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
                                            return ppipzsppideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                        } else {
                                            return 0;
                                        }
                                    case 1:
                                        if (c.getk() == 1) {
                                            if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(pzppi|pzppi)
                                                return ppipzppipzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                            } else {
                                                return 0;
                                            }
                                        } else {
                                            if (c.geti() == b.geti() && c.getj() == b.getj() && c.getk() == 0) {//(pzppi|ppi?)
                                                switch (d.getL()) {
                                                    case 0:
                                                        return ppipzsppideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                    case 1:
                                                        if (d.getk() == 1) {
                                                            return ppipzppipzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
                                        return sppisppideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                    } else {
                                        return 0;
                                    }
                                case 1:
                                    if (c.getk() == 1) {
                                        if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppis|pzppi)
                                            return sppippipzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                        } else {
                                            return 0;
                                        }
                                    } else {
                                        if (c.geti() == a.geti() && c.getj() == a.getj() && c.getk() == 0) {//(ppis|ppi?)
                                            switch (d.getL()) {
                                                case 0:
                                                    return sppisppideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                case 1:
                                                    if (d.getk() == 1) {
                                                        return sppippipzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
                                            return ppipzsppideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                        } else {
                                            return 0;
                                        }
                                    case 1:
                                        if (c.getk() == 1) {
                                            if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppipz|pzppi)
                                                return ppipzppipzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                            } else {
                                                return 0;
                                            }
                                        } else {
                                            if (c.geti() == a.geti() && c.getj() == a.getj() && c.getk() == 0) {//(ppipz|ppi?)
                                                switch (d.getL()) {
                                                    case 0:
                                                        return ppipzsppideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                    case 1:
                                                        if (d.getk() == 1) {
                                                            return ppipzppipzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
                                                    return ppippissderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                } else {
                                                    return 0;
                                                }
                                            case 1:
                                                if (d.getk() == 1 && a.geti() == b.geti() && a.getj() == b.getj()) {
                                                    return ppippispzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                } else {
                                                    return 0;
                                                }
                                        }

                                    case 1:
                                        if (c.getk() == 1) {
                                            switch (d.getL()) {
                                                case 0://(ppippi|pzs)
                                                    if (a.geti() == b.geti() && a.getj() == b.getj()) {
                                                        return ppippispzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                    } else {
                                                        return 0;
                                                    }

                                                case 1:
                                                    if (d.getk() == 1 && a.geti() == b.geti() && a.getj() == b.getj()) {//(ppippi|pzpz)
                                                        return ppippipzpzderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                    } else {
                                                        return 0;
                                                    }
                                            }
                                        } else {
                                            if (a.geti() == b.geti() && a.getj() == b.getj()) {//(pxpx|??) or (pypy|??)

                                                if (c.getL() == d.getL() && c.geti() == d.geti() && c.getj() == d.getj() && c.getk() == 0) {
                                                    if (a.geti() == c.geti()) {
                                                        return ppippippippideriv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                    } else {
                                                        return pxpxpypyderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
                                                    }
                                                } else {
                                                    return 0;
                                                }

                                            } else {//(pxpy|??) or (pypx|??)
                                                if (c.getL() == d.getL() && c.geti() != d.geti() && c.getj() != d.getj() && c.getk() == 0) {
                                                    return pxpypxpyderiv(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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

    public static double getGderiv(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau) {


        NDDO6G[] A = a.orbitalArray();
        NDDO6G[] B = b.orbitalArray();
        NDDO6G[] C = c.orbitalArray();
        NDDO6G[] D = d.orbitalArray();

        if (Math.abs(a.getCoords()[0] - c.getCoords()[0]) < 1E-3 && Math.abs(a.getCoords()[1] - c.getCoords()[1]) < 1E-3) {

            double[] coeffA2 = a.decomposition2(a.getCoords(), c.getCoords());
            double[] coeffB2 = b.decomposition2(a.getCoords(), c.getCoords());
            double[] coeffC2 = c.decomposition2(a.getCoords(), c.getCoords());
            double[] coeffD2 = d.decomposition2(a.getCoords(), c.getCoords());

            double[] coeffAderiv2 = derivativeDecomposition2(a.getCoords(), c.getCoords(), a, tau);
            double[] coeffBderiv2 = derivativeDecomposition2(a.getCoords(), c.getCoords(), b, tau);
            double[] coeffCderiv2 = derivativeDecomposition2(a.getCoords(), c.getCoords(), c, tau);
            double[] coeffDderiv2 = derivativeDecomposition2(a.getCoords(), c.getCoords(), d, tau);



            double sum2 = 0;

            for (int i = 0; i < coeffA2.length; i++) {
                for (int j = 0; j < coeffB2.length; j++) {
                    for (int k = 0; k < coeffC2.length; k++) {
                        for (int l = 0; l < coeffD2.length; l++) {


                            if (coeffA2[i] * coeffB2[j] * coeffC2[k] * coeffD2[l] != 0) {
                                sum2 += coeffA2[i] * coeffB2[j] * coeffC2[k] * coeffD2[l] * LocalTwoCenterERIderiv(A[i], B[j], C[k], D[l], tau) * 27.21;
                            }
                            if (coeffAderiv2[i] * coeffB2[j] * coeffC2[k] * coeffD2[l] != 0) {
                                sum2 += coeffAderiv2[i] * coeffB2[j] * coeffC2[k] * coeffD2[l] * NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * 27.21;
                            }

                            if (coeffA2[i] * coeffBderiv2[j] * coeffC2[k] * coeffD2[l] != 0) {
                                sum2 += coeffA2[i] * coeffBderiv2[j] * coeffC2[k] * coeffD2[l] * NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * 27.21;
                            }

                            if (coeffA2[i] * coeffB2[j] * coeffCderiv2[k] * coeffD2[l] != 0) {
                                sum2 += coeffA2[i] * coeffB2[j] * coeffCderiv2[k] * coeffD2[l] * NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * 27.21;
                            }

                            if (coeffA2[i] * coeffB2[j] * coeffC2[k] * coeffDderiv2[l] != 0) {
                                sum2 += coeffA2[i] * coeffB2[j] * coeffC2[k] * coeffDderiv2[l] * NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * 27.21;
                            }


                        }
                    }
                }
            }

            return sum2;
        }

        else {
            double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
            double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
            double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
            double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());

            double[] coeffAderiv = derivativeDecomposition(a.getCoords(), c.getCoords(), a, tau);
            double[] coeffBderiv = derivativeDecomposition(a.getCoords(), c.getCoords(), b, tau);
            double[] coeffCderiv = derivativeDecomposition(a.getCoords(), c.getCoords(), c, tau);
            double[] coeffDderiv = derivativeDecomposition(a.getCoords(), c.getCoords(), d, tau);



            double sum = 0;

            for (int i = 0; i < coeffA.length; i++) {
                for (int j = 0; j < coeffB.length; j++) {
                    for (int k = 0; k < coeffC.length; k++) {
                        for (int l = 0; l < coeffD.length; l++) {


                            if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
                                sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] * LocalTwoCenterERIderiv(A[i], B[j], C[k], D[l], tau) * 27.21;
                            }
                            if (coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
                                sum += coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] * NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * 27.21;
                            }

                            if (coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] != 0) {
                                sum += coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] * NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * 27.21;
                            }

                            if (coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] != 0) {
                                sum += coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] * NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * 27.21;
                            }

                            if (coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] != 0) {
                                sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] * NDDO6G.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * 27.21;
                            }


                        }
                    }
                }
            }




            return sum;
        }
    }

    public static double getGderivfinite(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau) {

        double orig = NDDO6G.getG(a, b, c, d);

        double[] newcoords = a.getCoords().clone();

        newcoords[tau] += 1E-8;

        NDDO6G anew = new NDDO6G(a, newcoords);
        NDDO6G bnew = new NDDO6G(b, newcoords);

        double perturbed = NDDO6G.getG(anew, bnew, c, d);

        return (perturbed - orig) / 1E-8;
    }



    public static double[] derivativeDecomposition(double[] point1, double[] point2, NDDO6G a, int tau) {

        if (a.getL() == 0) {
            return new double[]{0};
        }

        double R = GTO.R(point1, point2);
        double Rxy = Math.sqrt((point2[1] - point1[1]) * (point2[1] - point1[1]) + (point2[0] - point1[0]) * (point2[0] - point1[0]));


        switch (tau) {
            case 0:
                if (a.geti() == 1) {
                    double x1 = (point2[2] - point1[2]) / (R * Rxy) - (point2[0] - point1[0]) * (point2[0] - point1[0]) * (point2[2] - point1[2]) / (Rxy * Rxy * Rxy * R) - (point2[0] - point1[0]) * (point2[0] - point1[0]) * (point2[2] - point1[2]) / (R * R * R * Rxy);
                    double x2 = -(point2[0] - point1[0]) * (point2[1] - point1[1]) / (Rxy * Rxy * Rxy);
                    double x3 = (point2[0] - point1[0]) * (point2[0] - point1[0]) / (R * R * R) - 1 / R;

                    return new double[]{x1, x2, x3};
                } else if (a.getj() == 1) {
                    double x1 = -(point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2]) / (Rxy * Rxy * Rxy * R) - (point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2]) / (R * R * R * Rxy);
                    double x2 = (point2[0] - point1[0]) * (point2[0] - point1[0]) / (Rxy * Rxy * Rxy) - 1 / Rxy;
                    double x3 = (point2[0] - point1[0]) * (point2[1] - point1[1]) / (R * R * R);

                    return new double[]{x1, x2, x3};
                } else if (a.getk() == 1) {
                    double x1 = (point2[0] - point1[0]) * Rxy / (R * R * R) - (point2[0] - point1[0]) / (R * Rxy);
                    double x2 = 0;
                    double x3 = (point2[0] - point1[0]) * (point2[2] - point1[2]) / (R * R * R);

                    return new double[]{x1, x2, x3};
                }
            case 1:
                if (a.geti() == 1) {
                    double x1 = -(point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2]) / (Rxy * Rxy * Rxy * R) - (point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2]) / (R * R * R * Rxy);
                    double x2 = -(point2[1] - point1[1]) * (point2[1] - point1[1]) / (Rxy * Rxy * Rxy) + 1 / Rxy;
                    double x3 = (point2[0] - point1[0]) * (point2[1] - point1[1]) / (R * R * R);

                    return new double[]{x1, x2, x3};
                } else if (a.getj() == 1) {
                    double x1 = (point2[2] - point1[2]) / (R * Rxy) - (point2[2] - point1[2]) * (point2[1] - point1[1]) * (point2[1] - point1[1]) / (Rxy * Rxy * Rxy * R) - (point2[2] - point1[2]) * (point2[1] - point1[1]) * (point2[1] - point1[1]) / (R * R * R * Rxy);
                    double x2 = (point2[0] - point1[0]) * (point2[1] - point1[1]) / (Rxy * Rxy * Rxy);
                    double x3 = (point2[1] - point1[1]) * (point2[1] - point1[1]) / (R * R * R) - 1 / R;

                    return new double[]{x1, x2, x3};
                } else if (a.getk() == 1) {
                    double x1 = (point2[1] - point1[1]) * Rxy / (R * R * R) - (point2[1] - point1[1]) / (R * Rxy);
                    double x2 = 0;
                    double x3 = (point2[2] - point1[2]) * (point2[1] - point1[1]) / (R * R * R);

                    return new double[]{x1, x2, x3};
                }
            case 2:
                if (a.geti() == 1) {
                    double x1 = (point2[0] - point1[0]) / (R * Rxy) - (point2[0] - point1[0]) * (point2[2] - point1[2]) * (point2[2] - point1[2]) / (R * R * R * Rxy);
                    double x2 = 0;
                    double x3 = (point2[2] - point1[2]) * (point2[0] - point1[0]) / (R * R * R);

                    return new double[]{x1, x2, x3};
                } else if (a.getj() == 1) {
                    double x1 = (point2[1] - point1[1]) / (R * Rxy) - (point2[1] - point1[1]) * (point2[2] - point1[2]) * (point2[2] - point1[2]) / (R * R * R * Rxy);
                    double x2 = 0;
                    double x3 = (point2[2] - point1[2]) * (point2[1] - point1[1]) / (R * R * R);

                    return new double[]{x1, x2, x3};
                } else if (a.getk() == 1) {
                    double x1 = (point2[2] - point1[2]) * Rxy / (R * R * R);
                    double x2 = 0;
                    double x3 = (point2[2] - point1[2]) * (point2[2] - point1[2]) / (R * R * R) - 1 / R;

                    return new double[]{x1, x2, x3};
                }
        }
        return null;
    }

    public static double[] derivativeDecomposition2 (double[] point1, double[] point2, NDDO6G a, int tau) {

        if (a.getL() == 0) {
            return new double[] {0};
        }

        double x = point2[0] - point1[0];

        double y = point2[1] - point1[1];

        double z = point2[2] - point1[2];

        double Rxz = Math.sqrt(x * x + z * z);

        double R = GTO.R(point1, point2);

        if (a.getL() == 1) {
            switch (tau) {
                case 0:
                    if (a.geti() == 1) {
                        double val1 = x * x * y / (R * R * R * Rxz) + x * x * y / (R * Rxz * Rxz * Rxz) - y / (R * Rxz);
                        double val2 = - x * z / (Rxz * Rxz * Rxz);
                        double val3 = x * x / (R * R * R) - 1 / R;
                        return new double[] {val1, val2, val3};
                    }
                    else if (a.getj() == 1) {
                        double val1 = x / (R * Rxz) - x * Rxz / (R * R * R);
                        double val3 = x * y / (R * R * R);
                        return new double[] {val1, 0, val3};
                    }
                    else if (a.getk() == 1) {
                        double val1 = x * y * z / (R * R * R * Rxz) + x * y * z / (R * Rxz * Rxz * Rxz);
                        double val2 = x * x / (Rxz * Rxz * Rxz) - 1 / Rxz;
                        double val3 = x * z / (R * R * R);
                        return new double[] {val1, val2, val3};
                    }
                case 1:
                    if (a.geti() == 1) {
                        double val1 = x * y * y / (R * R * R * Rxz) - x / (R * Rxz);
                        double val3 = x * y / (R * R * R);
                        return new double[] {val1, 0, val3};
                    }
                    else if (a.getj() == 1) {
                        double val1 = - y * Rxz / (R * R * R);
                        double val3 = y * y / (R * R * R) - 1 / R;
                        return new double[] {val1, 0, val3};
                    }
                    else if (a.getk() == 1) {
                        double val1 = y * y * z / (R * R * R * Rxz) - z / (R * Rxz);
                        double val3 = y * z / (R * R * R);
                        return new double[] {val1, 0, val3};
                    }
                case 2:
                    if (a.geti() == 1) {
                        double val1 = x * y * z / (R * R * R * Rxz) + x * y * z / (R * Rxz * Rxz * Rxz);
                        double val2 = 1 / Rxz - z * z / (Rxz * Rxz * Rxz);
                        double val3 = x * z / (R * R * R);
                        return new double[] {val1, val2, val3};
                    }
                    else if (a.getj() == 1) {
                        double val1 = z / (R * Rxz) - z * Rxz / (R * R * R);
                        double val3 = z * y / (R * R * R);
                        return new double[] {val1, 0, val3};
                    }
                    else if (a.getk() == 1) {
                        double val1 = y * z * z / (R * R * R * Rxz) + y * z * z / (R * Rxz * Rxz * Rxz) - y / (R * Rxz);
                        double val2 = x * z / (Rxz * Rxz * Rxz);
                        double val3 = z * z / (R * R * R) - 1 / R;
                        return new double[] {val1, val2, val3};
                    }
            }
        }

        return null;

    }

    public static double gradient(NDDOAtom[] atoms, NDDOSolutionRestricted soln, int atomnum, int tau) {

        DoubleMatrix densitymatrix = soln.densityMatrix();


        NDDO6G[] orbitals = soln.orbitals;

        int[][] index = soln.index;

        int[] atomnumber = soln.atomNumber;

        DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;
                if (atomnumber[j] == atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int a = 0; a < atoms.length; a++) {
                            if (a != atomnum) {


                                sum -= atoms[a].getAtomProperties().getQ() * NDDODerivative.getGderiv(orbitals[j], orbitals[k], atoms[a].s(), atoms[a].s(), tau);
                            }
                        }
                    } else {
                        sum -= atoms[atomnum].getAtomProperties().getQ() * NDDODerivative.getGderiv(atoms[atomnum].s(), atoms[atomnum].s(), orbitals[j], orbitals[k], tau);
                    }
                } else {
                    if (atomnumber[j] == atomnum) {
                        sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) * LCGTO.getSDeriv(orbitals[j], orbitals[k], tau);
                    } else if (atomnumber[k] == atomnum) {
                        sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) * LCGTO.getSDeriv(orbitals[k], orbitals[j], tau);
                    }
                }

                H.put(j, k, sum);
                H.put(k, j, sum);
            }
        }

        DoubleMatrix G = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;
                if (atomnumber[j] == atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int a = 0; a < atoms.length; a++) {
                            if (a != atomnum) {

                                for (int l : index[a]) {
                                    if (l > -1) {
                                        for (int m : index[a]) {
                                            if (m > -1) {
                                                sum += densitymatrix.get(l, m) * NDDODerivative.getGderiv(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnum]) {
                                    if (m > -1) {
                                        sum += densitymatrix.get(l, m) * NDDODerivative.getGderiv(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);

                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (atomnumber[j] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[k]]) {
                                    if (m > -1) {
                                        sum -= 0.5 * densitymatrix.get(l, m) * NDDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    } else if (atomnumber[k] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[j]]) {
                                    if (m > -1) {
                                        sum -= 0.5 * densitymatrix.get(l, m) * NDDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    }
                }

                G.put(j, k, sum);
                G.put(k, j, sum);

            }
        }

        DoubleMatrix F = H.dup().add(G);

        double e = 0;

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = 0; k < orbitals.length; k++) {
                e += 0.5 * densitymatrix.get(j, k) * (H.get(j, k) + F.get(j, k));
            }
        }

        //System.out.println ("Electronic gradient (Analytic): " + e);

        for (int j = 0; j < atoms.length; j++) {
            if (j != atomnum) {
                e += atoms[atomnum].crfDeriv(atoms[j], tau);
            }
        }

        if (Math.abs(e - grad(atoms, soln, atomnum, tau)) > 1E-5) {
            System.err.println ("oh well, you can't say you weren't expecting that...");
            System.exit(0);
        }


        return e;

    }

    public static double gradientUnrestricted(NDDOAtom[] atoms, NDDOSolutionUnrestricted soln, int atomnum, int tau) {

        DoubleMatrix alphadensity = soln.alphaDensity();

        DoubleMatrix betadensity = soln.betaDensity();

        NDDO6G[] orbitals = soln.orbitals;

        int[][] index = soln.index;

        int[] atomnumber = soln.atomNumber;

        DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;
                if (atomnumber[j] == atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int a = 0; a < atoms.length; a++) {
                            if (a != atomnum) {


                                sum -= atoms[a].getAtomProperties().getQ() * NDDODerivative.getGderiv(orbitals[j], orbitals[k], atoms[a].s(), atoms[a].s(), tau);
                            }
                        }
                    } else {
                        sum -= atoms[atomnum].getAtomProperties().getQ() * NDDODerivative.getGderiv(atoms[atomnum].s(), atoms[atomnum].s(), orbitals[j], orbitals[k], tau);
                    }
                } else {
                    if (atomnumber[j] == atomnum) {
                        sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) * LCGTO.getSDeriv(orbitals[j], orbitals[k], tau);
                    } else if (atomnumber[k] == atomnum) {
                        sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) * LCGTO.getSDeriv(orbitals[k], orbitals[j], tau);
                    }
                }

                H.put(j, k, sum);
                H.put(k, j, sum);
            }
        }

        DoubleMatrix J = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;
                if (atomnumber[j] == atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int a = 0; a < atoms.length; a++) {
                            if (a != atomnum) {

                                for (int l : index[a]) {
                                    if (l > -1) {
                                        for (int m : index[a]) {
                                            if (m > -1) {
                                                sum += (alphadensity.get(l, m) + betadensity.get(l, m)) * NDDODerivative.getGderiv(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnum]) {
                                    if (m > -1) {
                                        sum += (alphadensity.get(l, m) + betadensity.get(l, m)) * NDDODerivative.getGderiv(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
                                    }
                                }
                            }
                        }
                    }
                }

                J.put(j, k, sum);
                J.put(k, j, sum);

            }
        }

        DoubleMatrix Ka = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;

                if (atomnumber[j] != atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[k]]) {
                                    if (m > -1) {
                                        sum -= alphadensity.get(l, m) * NDDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    } else if (atomnumber[k] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[j]]) {
                                    if (m > -1) {
                                        sum -= alphadensity.get(l, m) * NDDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    }
                }

                Ka.put(j, k, sum);
                Ka.put(k, j, sum);

            }
        }

        DoubleMatrix Kb = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;

                if (atomnumber[j] != atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[k]]) {
                                    if (m > -1) {
                                        sum -= betadensity.get(l, m) * NDDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    } else if (atomnumber[k] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[j]]) {
                                    if (m > -1) {
                                        sum -= betadensity.get(l, m) * NDDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    }
                }

                Kb.put(j, k, sum);
                Kb.put(k, j, sum);

            }
        }

        DoubleMatrix Fa = H.add(J).add(Ka);
        DoubleMatrix Fb = H.add(J).add(Kb);

        double e = 0;

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = 0; k < orbitals.length; k++) {
                e += 0.5 * alphadensity.get(j, k) * (H.get(j, k) + Fa.get(j, k));
                e += 0.5 * betadensity.get(j, k) * (H.get(j, k) + Fb.get(j, k));
            }
        }


        for (int j = 0; j < atoms.length; j++) {
            if (j != atomnum) {
                e += atoms[atomnum].crfDeriv(atoms[j], tau);
            }
        }

        if (Math.abs(e - grad(atoms, soln, atomnum, tau)) > 1E-5) {
            System.err.println ("oh well, you can't say you weren't expecting that...");
            System.exit(0);
        }

        return e;
    }

    public static DoubleMatrix[][] gradientroutine (NDDOAtom[] atoms, NDDOSolutionRestricted soln) {

        DoubleMatrix[] fockderivatives = new DoubleMatrix[atoms.length * 3];

        DoubleMatrix grad = new DoubleMatrix(atoms.length * 3, 1);

        int count = 0;

        for (int a = 0; a < atoms.length; a++) {
            for (int tau = 0; tau < 3; tau++) {
                DoubleMatrix[] matrices = staticderivs(atoms, soln, a, tau);
                fockderivatives[count] = matrices[1];
                double sum = 0;

                for (int i = 0; i < matrices[1].rows; i++) {
                    for (int j = 0; j < matrices[1].rows; j++) {
                        sum += 0.5 * soln.densityMatrix().get(i, j) * (matrices[0].get(i, j) + matrices[1].get(i, j));
                    }
                }

                for (int j = 0; j < atoms.length; j++) {
                    if (j != a) {
                        sum += atoms[a].crfDeriv(atoms[j], tau);
                    }
                }

                grad.put(count, 0, sum);

                count++;

            }
        }

        return new DoubleMatrix[][] {new DoubleMatrix[] {grad}, fockderivatives};
    }

    public static DoubleMatrix[][] gradientroutine (NDDOAtom[] atoms, NDDOSolutionUnrestricted soln) {

        DoubleMatrix[] fockderivativesalpha = new DoubleMatrix[atoms.length * 3];

        DoubleMatrix[] fockderivativesbeta = new DoubleMatrix[atoms.length * 3];

        DoubleMatrix grad = new DoubleMatrix(atoms.length * 3, 1);

        int count = 0;

        for (int a = 0; a < atoms.length; a++) {
            for (int tau = 0; tau < 3; tau++) {
                DoubleMatrix[] matrices = staticderivs(atoms, soln, a, tau);
                fockderivativesalpha[count] = matrices[1];
                fockderivativesbeta[count] = matrices[2];
                double sum = 0;

                for (int i = 0; i < matrices[1].rows; i++) {
                    for (int j = 0; j < matrices[1].rows; j++) {
                        sum += 0.5 * soln.densityMatrix().get(i, j) * (matrices[0].get(i, j) + matrices[1].get(i, j));
                        sum += 0.5 * soln.densityMatrix().get(i, j) * (matrices[0].get(i, j) + matrices[2].get(i, j));
                    }
                }

                for (int j = 0; j < atoms.length; j++) {
                    if (j != a) {
                        sum += atoms[a].crfDeriv(atoms[j], tau);
                    }
                }

                grad.put(count, 0, sum);

                count++;

            }
        }

        return new DoubleMatrix[][] {new DoubleMatrix[] {grad}, fockderivativesalpha, fockderivativesbeta};
    }

    public static double grad (NDDOAtom[] atoms, NDDOSolutionRestricted soln, int atomnum, int tau) {

        double e = 0;

        for (int a = 0; a < atoms.length; a++) {
            if (a != atomnum) {
                e += Ederiv (atomnum, a, soln.index, soln.densityMatrix(), atoms, soln.orbitals, tau);
                e += atoms[atomnum].crfDeriv(atoms[a], tau);
            }
        }

        return e;

    }

    public static double grad (NDDOAtom[] atoms, NDDOSolutionUnrestricted soln, int atomnum, int tau) {

        double e = 0;

        for (int a = 0; a < atoms.length; a++) {
            if (a != atomnum) {
                e += Ederiv (atomnum, a, soln.index, soln.alphaDensity(), soln.betaDensity(), atoms, soln.orbitals, tau);
                e += atoms[atomnum].crfDeriv(atoms[a], tau);
            }
        }

        return e;

    }

    private static double Ederiv( int atomnum1, int atomnum2, int[][] index, DoubleMatrix densityMatrix, NDDOAtom[] atoms, NDDO6G[] orbitals, int tau) {

        double e = 0;

        for (int i: index[atomnum1]) {
            for (int j: index[atomnum1]) {
                if (i != -1 && j != -1) {
                    e += densityMatrix.get(i, j) * atoms[atomnum2].Vderiv(orbitals[i], orbitals[j], tau);
                }
            }
        }

        for (int k: index[atomnum2]) {
            for (int l: index[atomnum2]) {
                if (k != -1 && l != -1) {
                    e -= densityMatrix.get(k, l) * atoms[atomnum1].Vderiv(orbitals[k], orbitals[l], tau);
                }
            }
        }

        for (int i: index[atomnum1]) {
            for (int k: index[atomnum2]) {
                if (i != -1 && k != -1) {
                    e += 2 * densityMatrix.get(i, k) * NDDO6G.betaderiv(orbitals[i], orbitals[k], tau);
                }
            }
        }

        for (int i: index[atomnum1]) {
            for (int j: index[atomnum1]) {
                for (int k: index[atomnum2]) {
                    for (int l: index[atomnum2]) {
                        if (i != -1 && j != -1 && k != -1 && l != -1) {
                            e +=(densityMatrix.get(i, j) * densityMatrix.get(k, l) - densityMatrix.get (i, k) * 0.5 * densityMatrix.get(j, l))
                                    * NDDODerivative.getGderiv(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau);
                        }
                    }
                }
            }
        }

        return e;


    }

    private static double Ederiv( int atomnum1, int atomnum2, int[][] index, DoubleMatrix alphaDensity, DoubleMatrix betaDensity, NDDOAtom[] atoms, NDDO6G[] orbitals, int tau) {

        double e = 0;


        for (int i: index[atomnum1]) {
            for (int j: index[atomnum1]) {
                if (i != -1 && j != -1) {
                    e += (alphaDensity.get(i, j) + betaDensity.get(i, j)) * atoms[atomnum2].Vderiv(orbitals[i], orbitals[j], tau);
                }
            }
        }

        for (int k: index[atomnum2]) {
            for (int l: index[atomnum2]) {
                if (k != -1 && l != -1) {
                    e -= (alphaDensity.get(k, l) + betaDensity.get(k, l)) * atoms[atomnum1].Vderiv(orbitals[k], orbitals[l], tau);
                }
            }
        }

        for (int i: index[atomnum1]) {
            for (int k: index[atomnum2]) {
                if (i != -1 && k != -1) {
                    e += 2 * (alphaDensity.get(i, k) + betaDensity.get(i, k)) * NDDO6G.betaderiv(orbitals[i], orbitals[k], tau);
                }
            }
        }

        for (int i: index[atomnum1]) {
            for (int j: index[atomnum1]) {
                for (int k: index[atomnum2]) {
                    for (int l: index[atomnum2]) {
                        if (i != -1 && j != -1 && k != -1 && l != -1) {
                            e +=((alphaDensity.get(i, j) + betaDensity.get(i, j)) * (alphaDensity.get(k, l) + betaDensity.get(k, l)) - alphaDensity.get (i, k) * alphaDensity.get(j, l) - betaDensity.get (i, k) * betaDensity.get(j, l))
                                    * NDDODerivative.getGderiv(orbitals[i], orbitals[j], orbitals[k], orbitals[l], tau);
                        }
                    }
                }
            }
        }

        return e;


    }

    public static DoubleMatrix[] staticderivs (NDDOAtom[] atoms, NDDOSolutionRestricted soln, int atomnum, int tau) {

        DoubleMatrix densitymatrix = soln.densityMatrix();


        NDDO6G[] orbitals = soln.orbitals;

        int[][] index = soln.index;

        int[] atomnumber = soln.atomNumber;

        DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;
                if (atomnumber[j] == atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int a = 0; a < atoms.length; a++) {
                            if (a != atomnum) {


                                sum -= atoms[a].getAtomProperties().getQ() * NDDODerivative.getGderiv(orbitals[j], orbitals[k], atoms[a].s(), atoms[a].s(), tau);
                            }
                        }
                    } else {
                        sum -= atoms[atomnum].getAtomProperties().getQ() * NDDODerivative.getGderiv(atoms[atomnum].s(), atoms[atomnum].s(), orbitals[j], orbitals[k], tau);
                    }
                } else {
                    if (atomnumber[j] == atomnum) {
                        sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) * LCGTO.getSDeriv(orbitals[j], orbitals[k], tau);
                    } else if (atomnumber[k] == atomnum) {
                        sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) * LCGTO.getSDeriv(orbitals[k], orbitals[j], tau);
                    }
                }

                H.put(j, k, sum);
                H.put(k, j, sum);
            }
        }

        DoubleMatrix G = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;
                if (atomnumber[j] == atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int a = 0; a < atoms.length; a++) {
                            if (a != atomnum) {

                                for (int l : index[a]) {
                                    if (l > -1) {
                                        for (int m : index[a]) {
                                            if (m > -1) {
                                                sum += densitymatrix.get(l, m) * NDDODerivative.getGderiv(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnum]) {
                                    if (m > -1) {
                                        sum += densitymatrix.get(l, m) * NDDODerivative.getGderiv(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);

                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (atomnumber[j] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[k]]) {
                                    if (m > -1) {
                                        sum -= 0.5 * densitymatrix.get(l, m) * NDDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    } else if (atomnumber[k] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[j]]) {
                                    if (m > -1) {
                                        sum -= 0.5 * densitymatrix.get(l, m) * NDDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    }
                }

                G.put(j, k, sum);
                G.put(k, j, sum);

            }
        }

        DoubleMatrix F = H.dup().add(G);

        return new DoubleMatrix[] {H, F};

    }

    public static DoubleMatrix[] staticderivs (NDDOAtom[] atoms, NDDOSolutionUnrestricted soln, int atomnum, int tau) {

        DoubleMatrix alphadensity = soln.alphaDensity();

        DoubleMatrix betadensity = soln.betaDensity();

        NDDO6G[] orbitals = soln.orbitals;

        int[][] index = soln.index;

        int[] atomnumber = soln.atomNumber;

        DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;
                if (atomnumber[j] == atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int a = 0; a < atoms.length; a++) {
                            if (a != atomnum) {


                                sum -= atoms[a].getAtomProperties().getQ() * NDDODerivative.getGderiv(orbitals[j], orbitals[k], atoms[a].s(), atoms[a].s(), tau);
                            }
                        }
                    } else {
                        sum -= atoms[atomnum].getAtomProperties().getQ() * NDDODerivative.getGderiv(atoms[atomnum].s(), atoms[atomnum].s(), orbitals[j], orbitals[k], tau);
                    }
                } else {
                    if (atomnumber[j] == atomnum) {
                        sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) * LCGTO.getSDeriv(orbitals[j], orbitals[k], tau);
                    } else if (atomnumber[k] == atomnum) {
                        sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) * LCGTO.getSDeriv(orbitals[k], orbitals[j], tau);
                    }
                }

                H.put(j, k, sum);
                H.put(k, j, sum);
            }
        }

        DoubleMatrix J = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;
                if (atomnumber[j] == atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int a = 0; a < atoms.length; a++) {
                            if (a != atomnum) {

                                for (int l : index[a]) {
                                    if (l > -1) {
                                        for (int m : index[a]) {
                                            if (m > -1) {
                                                sum += (alphadensity.get(l, m) + betadensity.get(l, m)) * NDDODerivative.getGderiv(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnum]) {
                                    if (m > -1) {
                                        sum += (alphadensity.get(l, m) + betadensity.get(l, m)) * NDDODerivative.getGderiv(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
                                    }
                                }
                            }
                        }
                    }
                }

                J.put(j, k, sum);
                J.put(k, j, sum);

            }
        }

        DoubleMatrix Ka = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;

                if (atomnumber[j] != atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[k]]) {
                                    if (m > -1) {
                                        sum -= alphadensity.get(l, m) * NDDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    } else if (atomnumber[k] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[j]]) {
                                    if (m > -1) {
                                        sum -= alphadensity.get(l, m) * NDDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    }
                }

                Ka.put(j, k, sum);
                Ka.put(k, j, sum);

            }
        }

        DoubleMatrix Kb = DoubleMatrix.zeros(orbitals.length, orbitals.length);

        for (int j = 0; j < orbitals.length; j++) {
            for (int k = j; k < orbitals.length; k++) {
                double sum = 0;

                if (atomnumber[j] != atomnumber[k]) {
                    if (atomnumber[j] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[k]]) {
                                    if (m > -1) {
                                        sum -= betadensity.get(l, m) * NDDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    } else if (atomnumber[k] == atomnum) {
                        for (int l : index[atomnum]) {
                            if (l > -1) {
                                for (int m : index[atomnumber[j]]) {
                                    if (m > -1) {
                                        sum -= betadensity.get(l, m) * NDDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
                                    }
                                }
                            }
                        }
                    }
                }

                Kb.put(j, k, sum);
                Kb.put(k, j, sum);

            }
        }

        DoubleMatrix Fa = H.add(J).add(Ka);
        DoubleMatrix Fb = H.add(J).add(Kb);

        return new DoubleMatrix[] {H, Fa, Fb};

    }

    public static DoubleMatrix densitymatrixderivfinite (NDDOAtom[] atoms,NDDOSolutionRestricted soln, int atomnum, int tau) {

        DoubleMatrix orig = soln.densityMatrix();



        NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum, tau);

        DoubleMatrix perturbed = new NDDOSolutionRestricted (newatoms, soln.charge).densityMatrix();

        return perturbed.sub(orig).mmul(1E7);


    }

    public static DoubleMatrix[] densitymatrixderivfinite (NDDOAtom[] atoms,NDDOSolutionUnrestricted soln, int atomnum, int tau) {

        DoubleMatrix aorig = soln.alphaDensity();

        DoubleMatrix borig = soln.betaDensity();

        NDDOAtom[] newatoms = Utils.perturbAtomCoords(atoms, atomnum, tau);

        NDDOSolutionUnrestricted newsoln = new NDDOSolutionUnrestricted (newatoms, soln.charge, soln.multiplicity);

        DoubleMatrix aperturbed = newsoln.alphaDensity();

        DoubleMatrix bperturbed = newsoln.betaDensity();

        return new DoubleMatrix[] {aperturbed.sub(aorig).mmul(1E7), bperturbed.sub(borig).mmul(1E7)};


    }


}
