package mndoparam.mndo;

import org.jblas.DoubleMatrix;

import scf.GTO;
import scf.LCGTO;

public class MNDODerivative {

	
	private static double qqderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a00 = p01 + p02;
		return (xB[tau] - xA[tau]) * Math.pow(R * R + a00 * a00, -1.5);
	}
	
	private static double quzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a01 = p01 + p12;
		return 0.5/R * (xB[tau] - xA[tau]) * (R + D12) * Math.pow((R + D12) * (R + D12) + a01 * a01, -1.5) 
			 - 0.5/R * (xB[tau] - xA[tau]) * (R - D12) * Math.pow((R - D12) * (R - D12) + a01 * a01, -1.5);
	}
	
	private static double qQpipideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a02 = p01 + p22;
		return 0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D22 * D22 + a02 * a02, -1.5)
			  -0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + a02 * a02, -1.5);
	}
	private static double qQzzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a02 = p01 + p22;
		return 0.25/R * (xB[tau] - xA[tau]) * (R + 2 * D22) * Math.pow((R + 2 * D22) * (R + 2 * D22) + a02 * a02, -1.5) 
			  +0.25/R * (xB[tau] - xA[tau]) * (R - 2 * D22) * Math.pow((R - 2 * D22) * (R - 2 * D22) + a02 * a02, -1.5) 
			  - 0.5 * (xB[tau] - xA[tau])  * Math.pow(R * R + a02 * a02, -1.5);
	}
	
	private static double upiupideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a11 = p11 + p12;
		return 0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + (D11-D12) * (D11-D12) + a11 * a11, -1.5)
			  -0.5 * (xB[tau] - xA[tau]) * Math.pow(R * R + (D11+D12) * (D11+D12) + a11 * a11, -1.5);
	}
	
	private static double uzuzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a11 = p11 + p12;
		return 0.25/R * (xB[tau] - xA[tau]) * (R + D11 - D12) * Math.pow((R + D11 - D12) * (R + D11 - D12) + a11 * a11, -1.5)
			  -0.25/R * (xB[tau] - xA[tau]) * (R + D11 + D12) * Math.pow((R + D11 + D12) * (R + D11 + D12) + a11 * a11, -1.5)
			  -0.25/R * (xB[tau] - xA[tau]) * (R - D11 - D12) * Math.pow((R - D11 - D12) * (R - D11 - D12) + a11 * a11, -1.5)
			  +0.25/R * (xB[tau] - xA[tau]) * (R - D11 + D12) * Math.pow((R - D11 + D12) * (R - D11 + D12) + a11 * a11, -1.5);
	}
	
	private static double upiQpizderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a12 = p11 + p22;
		return -0.25/R * (xB[tau] - xA[tau]) * (R - D22) * Math.pow((R - D22) * (R - D22) + (D11-D22) * (D11-D22) + a12 * a12, -1.5)
			   +0.25/R * (xB[tau] - xA[tau]) * (R - D22) * Math.pow((R - D22) * (R - D22) + (D11+D22) * (D11+D22) + a12 * a12, -1.5)
			   +0.25/R * (xB[tau] - xA[tau]) * (R + D22) * Math.pow((R + D22) * (R + D22) + (D11-D22) * (D11-D22) + a12 * a12, -1.5)
			   -0.25/R * (xB[tau] - xA[tau]) * (R + D22) * Math.pow((R + D22) * (R + D22) + (D11+D22) * (D11+D22) + a12 * a12, -1.5);
	}
	
	private static double uzQpipideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a12 = p11 + p22;
		return -0.25/R * (xB[tau] - xA[tau]) * (R + D11) * Math.pow((R + D11) * (R + D11) + 4 * D22 * D22 + a12 * a12, -1.5)
			   +0.25/R * (xB[tau] - xA[tau]) * (R - D11) * Math.pow((R - D11) * (R - D11) + 4 * D22 * D22 + a12 * a12, -1.5)
			   +0.25/R * (xB[tau] - xA[tau]) * (R + D11) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
			   -0.25/R * (xB[tau] - xA[tau]) * (R - D11) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5);
	}
	
	private static double uzQzzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a12 = p11 + p22;
		return -0.125/R * (xB[tau] - xA[tau]) * (R + D11 - 2 * D22) * Math.pow((R + D11 - 2 * D22) * (R + D11 - 2 * D22) + a12 * a12, -1.5)
			   +0.125/R * (xB[tau] - xA[tau]) * (R - D11 - 2 * D22) * Math.pow((R - D11 - 2 * D22) * (R - D11 - 2 * D22) + a12 * a12, -1.5)
			   -0.125/R * (xB[tau] - xA[tau]) * (R + D11 + 2 * D22) * Math.pow((R + D11 + 2 * D22) * (R + D11 + 2 * D22) + a12 * a12, -1.5)
			   +0.125/R * (xB[tau] - xA[tau]) * (R - D11 + 2 * D22) * Math.pow((R - D11 + 2 * D22) * (R - D11 + 2 * D22) + a12 * a12, -1.5)
			   + 0.25/R * (xB[tau] - xA[tau]) * (R + D11) * Math.pow((R + D11) * (R + D11) + a12 * a12, -1.5)
			   - 0.25/R * (xB[tau] - xA[tau]) * (R - D11) * Math.pow((R - D11) * (R - D11) + a12 * a12, -1.5); 
	}
	
	private static double QpipiQpipideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.125 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
			  +0.125 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * (D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
			  - 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
			  - 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5)
			  + 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
	}
	
	private static double QxxQyyderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + 4 * D22 * D22 + a22 * a22, -1.5)
			  -0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
			  -0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D22 * D22 + a22 * a22, -1.5)
			  +0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
	}
	
	private static double QpipiQzzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.125/R * (xB[tau] - xA[tau]) * (R - 2 * D22) * Math.pow((R - 2 * D22) * (R - 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5)
			  +0.125/R * (xB[tau] - xA[tau]) * (R + 2 * D22) * Math.pow((R + 2 * D22) * (R + 2 * D22) + 4 * D21 * D21 + a22 * a22, -1.5)
			  -0.125/R * (xB[tau] - xA[tau]) * (R - 2 * D22) * Math.pow((R - 2 * D22) * (R - 2 * D22) + a22 * a22, -1.5)
			  -0.125/R * (xB[tau] - xA[tau]) * (R + 2 * D22) * Math.pow((R + 2 * D22) * (R + 2 * D22) + a22 * a22, -1.5)
			  -0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + 4 * D21 * D21 + a22 * a22, -1.5)
			  +0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
	}
	
	private static double QzzQzzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.0625/R * (xB[tau] - xA[tau]) * (R + 2*D21 - 2*D22) * Math.pow((R + 2*D21 - 2*D22) * (R + 2*D21 - 2*D22) + a22 * a22, -1.5)
			  +0.0625/R * (xB[tau] - xA[tau]) * (R + 2*D21 + 2*D22) * Math.pow((R + 2*D21 + 2*D22) * (R + 2*D21 + 2*D22) + a22 * a22, -1.5)
			  +0.0625/R * (xB[tau] - xA[tau]) * (R - 2*D21 - 2*D22) * Math.pow((R - 2*D21 - 2*D22) * (R - 2*D21 - 2*D22) + a22 * a22, -1.5)
			  +0.0625/R * (xB[tau] - xA[tau]) * (R - 2*D21 + 2*D22) * Math.pow((R - 2*D21 + 2*D22) * (R - 2*D21 + 2*D22) + a22 * a22, -1.5)
			  - 0.125/R * (xB[tau] - xA[tau]) * (R + 2*D21) * Math.pow((R + 2*D21) * (R + 2*D21) + a22 * a22, -1.5)
			  - 0.125/R * (xB[tau] - xA[tau]) * (R - 2*D21) * Math.pow((R - 2*D21) * (R - 2*D21) + a22 * a22, -1.5)
			  - 0.125/R * (xB[tau] - xA[tau]) * (R + 2*D22) * Math.pow((R + 2*D22) * (R + 2*D22) + a22 * a22, -1.5)
			  - 0.125/R * (xB[tau] - xA[tau]) * (R - 2*D22) * Math.pow((R - 2*D22) * (R - 2*D22) + a22 * a22, -1.5)
			  +  0.25 * (xB[tau] - xA[tau]) * Math.pow(R * R + a22 * a22, -1.5);
	}
	
	private static double QpizQpizderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		
		double R = GTO.R(xA, xB);
		double a22 = p21 + p22;
		return 0.125/R * (xB[tau] - xA[tau]) * (R + D21 - D22) * Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
			  -0.125/R * (xB[tau] - xA[tau]) * (R + D21 - D22) * Math.pow((R + D21 - D22) * (R + D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
			  -0.125/R * (xB[tau] - xA[tau]) * (R + D21 + D22) * Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
			  +0.125/R * (xB[tau] - xA[tau]) * (R + D21 + D22) * Math.pow((R + D21 + D22) * (R + D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
			  -0.125/R * (xB[tau] - xA[tau]) * (R - D21 - D22) * Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
			  +0.125/R * (xB[tau] - xA[tau]) * (R - D21 - D22) * Math.pow((R - D21 - D22) * (R - D21 - D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5)
			  +0.125/R * (xB[tau] - xA[tau]) * (R - D21 + D22) * Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 - D22) * (D21 - D22) + a22 * a22, -1.5)
			  -0.125/R * (xB[tau] - xA[tau]) * (R - D21 + D22) * Math.pow((R - D21 + D22) * (R - D21 + D22) + (D21 + D22) * (D21 + D22) + a22 * a22, -1.5);
	}
	private static double ssssderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double ssppippideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double sspzpzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double ppippissderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau); 
	}
	
	private static double pzpzssderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}
	private static double ppippippippideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QpipiQpipideriv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double pxpxpypyderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QxxQyyderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double ppippipzpzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QpipiQzzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double pzpzppippideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQpipideriv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QpipiQzzderiv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}
	
	private static double pzpzpzpzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return qqderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) + qQzzderiv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + QzzQzzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double spzssderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return -quzderiv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}
	
	private static double spzppippideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return -quzderiv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + uzQpipideriv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double spzpzpzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return -quzderiv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau) + uzQzzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double ssspzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return quzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double ppippispzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return quzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) -uzQpipideriv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}
	
	private static double pzpzspzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return quzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) - uzQzzderiv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}
	
	private static double sppisppideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return upiupideriv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double spzspzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return uzuzderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double sppippipzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return upiQpizderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double ppipzsppideriv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return -upiQpizderiv (p02, p12, p22, D12, D22, p01, p11, p21, D11, D21, xA, xB, tau);
	}
	
	private static double ppipzppipzderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return QpizQpizderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau);
	}
	
	private static double pxpypxpyderiv (double p01, double p11, double p21, double D11, double D21, double p02, double p12, double p22, double D12, double D22, double[] xA, double[] xB, int tau) {
		return 0.5 * (ppippippippideriv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau) - pxpxpypyderiv (p01, p11, p21, D11, D21, p02, p12, p22, D12, D22, xA, xB, tau));
	}
	
	private static double LocalTwoCenterERIderiv (MNDO6G a, MNDO6G b, MNDO6G c, MNDO6G d, int tau) {
		
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
						return ssssderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
						
					case 1:
						if (d.getk() == 1) {//(ss|spz)
							return ssspzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
						}
						else {//(ss|sppi) = 0
							return 0;
						}
					default:
						System.err.println ("oh no");
						return 0;
					}
				
				case 1: //(ss|p?)
					if (c.getk() == 1) {//(ss|pz?)
						
						switch (d.getL()) {
						
						case 0://(ss|pzs)
							return ssspzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							
						case 1:
							if (d.getk() == 1) {//(ss|pzpz)
								return sspzpzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							}
							else {//(ss|pzppi) = 0
								return 0;
							}
						default:
							return 0;
						}
					}
					else {//(ss|ppi?)
						
						if (d.getL() == 1 && d.getk() == 0 && c.geti() == d.geti() && c.getj() == d.getj()) {//(ss|ppippi)
							return ssppippideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
						}
						else {//all others are 0
							return 0;
						}
					}
				default:
					System.err.println ("oh no");
					return 0;
					
				}
			case 1: //(sp|??)
				
				if (b.getk() == 1) {//(spz|??)
					
					switch (c.getL()) {
					
					case 0://(spz|s?)
						
						switch (d.getL()) {
						
						case 0://(spz|ss)
							return spzssderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							
						case 1:
							if (d.getk() == 1) {//(spz|spz)
								return spzspzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							}
							else {
								return 0;
							}
						}
					
					case 1:
						if (c.getk() == 1) {//(spz|pz?)
							
							switch (d.getL()) {
							
							case 0://(spz|pzs)
								return spzspzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							
							case 1:
								if (d.getk() == 1) {//(spz|pzpz)
									return spzpzpzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								}
								else {//(spz|pzppi) = 0
									return 0;
								}
							}
						}
						else {//(spz|ppi?)
							if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
								return spzppippideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							}
							else {
								return 0;
							}
						}
					default:
						System.err.println ("oh no");
						return 0;
					}
				}
				else {//(sppi|??)
					
					switch (c.getL()) {
					case 0://(sppi|s?)
						if (d.geti()==b.geti() && d.getj()==b.getj() && d.getk() == 0) {//(sppi|sppi)
							return sppisppideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
						}
						else {
							return 0;
						}
					case 1:
						if (c.getk() == 1) {
							if (d.geti()==b.geti() && d.getj()==b.getj() && d.getk() == 0) {//(sppi|pzppi)
								return sppippipzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							}
							else {
								return 0;
							}
						}
						else {
							if (c.geti()==b.geti() && c.getj()==b.getj() && c.getk() == 0) {//(sppi|ppi?)
								switch (d.getL()) {
								case 0:
									return sppisppideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								case 1:
									if (d.getk() == 1) {
										return sppippipzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
									}
									else {
										return 0;
									}
								default:
									return 0;
								}
							}
							else {
								return 0;
							}
						}
					default:
						System.err.println ("oh no");
						return 0;
					}
				}
			default:
				System.err.println ("oh no");
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
							return spzssderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							
						case 1:
							if (d.getk() == 1) {//(pzs|spz)
								return spzspzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							}
							else {
								return 0;
							}
						}
					
					case 1:
						if (c.getk() == 1) {//(pzs|pz?)
							
							switch (d.getL()) {
							
							case 0://(pzs|pzs)
								return spzspzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							
							case 1:
								if (d.getk() == 1) {//(pzs|pzpz)
									return spzpzpzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								}
								else {//(pzs|pzppi) = 0
									return 0;
								}
							}
						}
						else {//(pzs|ppi?)
							if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
								return spzppippideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							}
							else {
								return 0;
							}
						}
					default:
						System.err.println ("oh no");
						return 0;
					}
				case 1:
					
					if (b.getk() == 1) {//(pzpz|??)
						
						switch (c.getL()) {
						
						case 0://(pzpz|s?)
							
							switch (d.getL()) {
							
							case 0://(pzpz|ss)
								return pzpzssderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								
							case 1:
								if (d.getk() == 1) {//(pzpz|spz)
									return pzpzspzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								}
								else {
									return 0;
								}
							}
						
						case 1:
							if (c.getk() == 1) {//(pzpz|pz?)
								
								switch (d.getL()) {
								
								case 0://(pzpz|pzs)
									return pzpzspzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								
								case 1:
									if (d.getk() == 1) {//(pzpz|pzpz)
										return pzpzpzpzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
									}
									else {//(pzpz|pzppi) = 0
										return 0;
									}
								}
							}
							else {//(pzpz|ppi?)
								if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
									return pzpzppippideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								}
								else {
									return 0;
								}
							}
						default:
							System.err.println ("oh no");
							return 0;
						}
					}
					else {//(pzppi|??)
						
						switch (c.getL()) {
						case 0://(pzppi|s?)
							if (d.geti()==b.geti() && d.getj()==b.getj() && d.getk() == 0) {//(pzppi|sppi)
								return ppipzsppideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							}
							else {
								return 0;
							}
						case 1:
							if (c.getk() == 1) {
								if (d.geti()==b.geti() && d.getj()==b.getj() && d.getk() == 0) {//(pzppi|pzppi)
									return ppipzppipzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								}
								else {
									return 0;
								}
							}
							else {
								if (c.geti()==b.geti() && c.getj()==b.getj() && c.getk() == 0) {//(pzppi|ppi?)
									switch (d.getL()) {
									case 0:
										return ppipzsppideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
									case 1:
										if (d.getk() == 1) {
											return ppipzppipzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
										}
										else {
											return 0;
										}
									default:
										return 0;
									}
								}
								else {
									return 0;
								}
							}
						default:
							System.err.println ("oh no");
							return 0;
						}
					}
				}
			}
			else {//(ppi?|??);
				
				switch (b.getL()) {
				case 0://(ppis|??)
					
					switch (c.getL()) {
					case 0://(ppis|s?)
						if (d.geti()==a.geti() && d.getj()==a.getj() && d.getk() == 0) {//(ppis|sppi)
							return sppisppideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
						}
						else {
							return 0;
						}
					case 1:
						if (c.getk() == 1) {
							if (d.geti()==a.geti() && d.getj()==a.getj() && d.getk() == 0) {//(ppis|pzppi)
								return sppippipzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							}
							else {
								return 0;
							}
						}
						else {
							if (c.geti()==a.geti() && c.getj()==a.getj() && c.getk() == 0) {//(ppis|ppi?)
								switch (d.getL()) {
								case 0:
									return sppisppideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								case 1:
									if (d.getk() == 1) {
										return sppippipzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
									}
									else {
										return 0;
									}
								default:
									return 0;
								}
							}
							else {
								return 0;
							}
						}
					default:
						System.err.println ("oh no");
						return 0;
					}
				case 1:
					if (b.getk() == 1) {//(ppipz|??)
						switch (c.getL()) {
						case 0://(ppipz|s?)
							if (d.geti()==a.geti() && d.getj()==a.getj() && d.getk() == 0) {//(ppipz|sppi)
								return ppipzsppideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
							}
							else {
								return 0;
							}
						case 1:
							if (c.getk() == 1) {
								if (d.geti()==a.geti() && d.getj()==a.getj() && d.getk() == 0) {//(ppipz|pzppi)
									return ppipzppipzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								}
								else {
									return 0;
								}
							}
							else {
								if (c.geti()==a.geti() && c.getj()==a.getj() && c.getk() == 0) {//(ppipz|ppi?)
									switch (d.getL()) {
									case 0:
										return ppipzsppideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
									case 1:
										if (d.getk() == 1) {
											return ppipzppipzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
										}
										else {
											return 0;
										}
									default:
										return 0;
									}
								}
								else {
									return 0;
								}
							}
						default:
							System.err.println ("oh no");
							return 0;
						}
						
					}
					else {
						
						switch (c.getL()) {
						case 0://(ppippi|s?)
							switch (d.getL()) {
							case 0://(ppippi|ss)
								if (a.geti() == b.geti() && a.getj() == b.getj()) {
									return ppippissderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								}
								else {
									return 0;
								}
							case 1:
								if (d.getk() ==1&&a.geti() == b.geti() && a.getj() == b.getj()) {
									return ppippispzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
								}
								else {
									return 0;
								}
							}
							
						case 1:
							if (c.getk() == 1) {
								switch (d.getL()) {
								case 0://(ppippi|pzs)
									if (a.geti() == b.geti() && a.getj() == b.getj()) {
										return ppippispzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
									}
									else {
										return 0;
									}
									
								case 1:
									if (d.getk() ==1&&a.geti() == b.geti() && a.getj() == b.getj()) {//(ppippi|pzpz)
										return ppippipzpzderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
									}
									else {
										return 0;
									}
								}
							}
							else {
								if (a.geti() == b.geti() && a.getj() == b.getj()) {//(pxpx|??) or (pypy|??)
									
									if (c.getL() == d.getL() && c.geti() == d.geti() && c.getj() == d.getj() && c.getk() == 0) {
										if (a.geti() == c.geti()) {
											return ppippippippideriv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
										}
										else {
											return pxpxpypyderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
										}
									}
									else {
										return 0;
									}
									
								}
								else {//(pxpy|??) or (pypx|??)
									if (c.getL() == d.getL() && c.geti() != d.geti() && c.getj() != d.getj() && c.getk() == 0) {
										return pxpypxpyderiv (a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, A, C, tau);
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
	
	public static double getGderiv (MNDO6G a, MNDO6G b, MNDO6G c, MNDO6G d, int tau) {
		
		if (a.getCoords()[0] - c.getCoords()[0] == 0 && a.getCoords()[1] - c.getCoords()[1] == 0) {
			return getGderivfinite (a,  b,  c,  d, tau);
		}
		
		double[] coeffA = a.decomposition(a.getCoords(), c.getCoords());
		double[] coeffB = b.decomposition(a.getCoords(), c.getCoords());
		double[] coeffC = c.decomposition(a.getCoords(), c.getCoords());
		double[] coeffD = d.decomposition(a.getCoords(), c.getCoords());
		
		double[] coeffAderiv = derivativedecomposition (a.getCoords(), c.getCoords(), a, tau);
		double[] coeffBderiv = derivativedecomposition (a.getCoords(), c.getCoords(), b, tau);
		double[] coeffCderiv = derivativedecomposition (a.getCoords(), c.getCoords(), c, tau);
		double[] coeffDderiv = derivativedecomposition (a.getCoords(), c.getCoords(), d, tau);
		
		
		MNDO6G[] A = a.orbitalarray2();
		MNDO6G[] B = b.orbitalarray2();
		MNDO6G[] C = c.orbitalarray2();
		MNDO6G[] D = d.orbitalarray2();
		
		double sum = 0;
		
		for (int i = 0; i < coeffA.length; i++) {
			for (int j = 0; j < coeffB.length; j++) {
				for (int k = 0; k < coeffC.length; k++) {
					for (int l = 0; l < coeffD.length; l++) {
						
						
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] * LocalTwoCenterERIderiv (A[i], B[j], C[k], D[l], tau) * 27.21;
						}
						if (coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffAderiv[i] * coeffB[j] * coeffC[k] * coeffD[l] * MNDO6G.LocalTwoCenterERI (A[i], B[j], C[k], D[l]) * 27.21;
						}
						
						if (coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffBderiv[j] * coeffC[k] * coeffD[l] * MNDO6G.LocalTwoCenterERI (A[i], B[j], C[k], D[l]) * 27.21;
						}
						
						if (coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffCderiv[k] * coeffD[l] * MNDO6G.LocalTwoCenterERI (A[i], B[j], C[k], D[l]) * 27.21;
						}
						
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] != 0) {
							sum += coeffA[i] * coeffB[j] * coeffC[k] * coeffDderiv[l] * MNDO6G.LocalTwoCenterERI (A[i], B[j], C[k], D[l]) * 27.21;
						}
						
						
					}
				}
			}
		}
		
		
		return sum;
	}
	
	public static double[] derivativedecompositionfinite (double[] point1, double[] point2, MNDO6G a, int tau) {
		if (a.getL()== 0) {
			return new double[] {0};
		}
		
		double[] orig = a.decomposition(point1, point2);
		
		point1 = point1.clone();
		
		point1[tau] += 1E-9;
		
		double[] perturbed = a.decomposition(point1, point2);
		
		return new double[] {(perturbed[0] - orig[0])/1E-9, (perturbed[1] - orig[1])/1E-9, (perturbed[2] - orig[2])/1E-9};
	}
	
	public static double[] derivativedecomposition (double[] point1, double[] point2, MNDO6G a, int tau) {
		
		if (a.getL() == 0) {
			return new double[] {0};
		}
		
		double R = GTO.R(point1, point2);
		double Rxy = Math.sqrt((point2[1] - point1[1]) * (point2[1] - point1[1]) + (point2[0] - point1[0]) * (point2[0] - point1[0]));
		
		
		switch (tau) {
		case 0:
			if (a.geti()== 1) {
				double x1 = (point2[2] - point1[2])/(R* Rxy) - (point2[0] - point1[0]) * (point2[0] - point1[0]) * (point2[2] - point1[2])/(Rxy * Rxy * Rxy * R) - (point2[0] - point1[0]) * (point2[0] - point1[0]) * (point2[2] - point1[2])/(R * R * R * Rxy);
				double x2 = -(point2[0] - point1[0]) * (point2[1] - point1[1])/(Rxy * Rxy * Rxy);
				double x3 =  (point2[0] - point1[0]) * (point2[0] - point1[0])/(R * R * R) - 1/R;
				
				return new double[] {x1, x2, x3};
			}
			else if (a.getj()==1) {
				double x1 = -(point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2])/(Rxy * Rxy * Rxy * R) - (point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2])/(R * R * R * Rxy);
				double x2 = (point2[0] - point1[0]) * (point2[0] - point1[0]) / (Rxy * Rxy * Rxy) - 1/Rxy;
				double x3 = (point2[0] - point1[0]) * (point2[1] - point1[1]) / (R * R * R);
				
				return new double[] {x1, x2, x3};
			}
			else if (a.getk() == 1) {
				double x1 = (point2[0] - point1[0]) * Rxy / (R * R * R) - (point2[0] - point1[0])/(R * Rxy);
				double x2 = 0;
				double x3 = (point2[0] - point1[0]) * (point2[2] - point1[2]) / (R * R * R);
				
				return new double[] {x1, x2, x3};
			}
		case 1:
			if (a.geti() == 1) {
				double x1 = -(point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2])/(Rxy * Rxy * Rxy * R) - (point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2])/(R * R * R * Rxy);
				double x2 = -(point2[1] - point1[1]) * (point2[1] - point1[1]) / (Rxy * Rxy * Rxy) + 1/Rxy;
				double x3 = (point2[0] - point1[0]) * (point2[1] - point1[1]) / (R * R * R);
				
				return new double[] {x1, x2, x3};
			}
			else if (a.getj() == 1) {
				double x1 = (point2[2] - point1[2])/(R* Rxy) - (point2[2] - point1[2]) * (point2[1] - point1[1]) * (point2[1] - point1[1])/(Rxy * Rxy * Rxy * R) - (point2[2] - point1[2]) * (point2[1] - point1[1]) * (point2[1] - point1[1])/(R * R * R * Rxy);
				double x2 = (point2[0] - point1[0]) * (point2[1] - point1[1])/(Rxy * Rxy * Rxy);
				double x3 = (point2[1] - point1[1]) * (point2[1] - point1[1])/(R * R * R) - 1/R;
				
				return new double[] {x1, x2, x3};
			}
			else if (a.getk() == 1) {
				double x1 = (point2[1] - point1[1]) * Rxy / (R * R * R) - (point2[1] - point1[1])/(R * Rxy);
				double x2 = 0;
				double x3 = (point2[2] - point1[2]) * (point2[1] - point1[1]) / (R * R * R);
				
				return new double[] {x1, x2, x3};
			}
		case 2:
			if (a.geti() == 1) {
				double x1 = (point2[0] - point1[0])/(R * Rxy) - (point2[0] - point1[0]) * (point2[2] - point1[2]) * (point2[2] - point1[2]) / (R * R * R * Rxy);
				double x2 = 0;
				double x3 = (point2[2] - point1[2]) * (point2[0] - point1[0]) / (R * R * R);
				
				return new double[] {x1, x2, x3};
			}
			else if (a.getj() == 1) {
				double x1 = (point2[1] - point1[1])/(R * Rxy) - (point2[1] - point1[1]) * (point2[2] - point1[2]) * (point2[2] - point1[2]) / (R * R * R * Rxy);
				double x2 = 0;
				double x3 = (point2[2] - point1[2]) * (point2[1] - point1[1]) / (R * R * R);
				
				return new double[] {x1, x2, x3};
			}
			else if (a.getk() == 1) {
				double x1 = (point2[2] - point1[2]) * Rxy / (R * R * R);
				double x2 = 0;
				double x3 = (point2[2] - point1[2]) * (point2[2] - point1[2])/(R * R * R) - 1/R;
				
				return new double[] {x1, x2, x3};
			}
		}
		return null;
	}

	
	public static double crfderivfinite (MNDOAtom a, MNDOAtom b, int tau) {
		return MNDOAtom.crfDeriv(a, b, tau);
	}
	
	public static double gradient (MNDOAtom[] atoms, DoubleMatrix densitymatrix, int atomnum, int tau) {
		//System.err.println ("Gradient evaluating");
		
		int i = 0;
		
		for (MNDOAtom a: atoms) {
			i += a.getOrbitals().length;
		}
		
		MNDO6G[] orbitals = new MNDO6G[i];
		
		i = 0;
		
		int[][] index = new int[atoms.length][4];
		int[] atomnumber = new int [orbitals.length];
		int count = 0;
		int count2 = 0;
		for (MNDOAtom a: atoms) {
			count2 = 0;
			for (MNDO6G orbital: a.getOrbitals()) {
				orbitals[i] = orbital;
				index[count][count2] = i;
				atomnumber[i] = count;
				i++;
				count2++;
			}
			
			
			if (a.getZ() == 1) {
				index[count][1] = -1;
				index[count][2] = -1;
				index[count][3] = -1;
				
			}
			count++;
		}
		
		int[][] missingindex = new int[atoms.length][4 * atoms.length - 4];
		
		for (int j = 0; j < atoms.length; j++) {
			for (int k = 0; k < 4 * atoms.length - 4; k++) {
				missingindex[j][k] = -1;
			}
		}
		
		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[] {index[j][0], index[j][1], index[j][2], index[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (nums[0] != k && nums[1] != k && nums[2] != k && nums[3] != k) {
					missingindex[j][counter] = k;
					counter++;
				}
			}
		}
		
		DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {
								
								
								sum -= atoms[a].getQ() * MNDODerivative.getGderiv(orbitals[j], orbitals[k], atoms[a].s(), atoms[a].s(), tau);
							}
						}
					}
					else {
						sum -= atoms[atomnum].getQ() * MNDODerivative.getGderiv(atoms[atomnum].s(), atoms[atomnum].s(), orbitals[j], orbitals[k], tau);
					}
				}
				
				else {
					if (atomnumber[j] == atomnum) {
						sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) *LCGTO.getSderiv(orbitals[j], orbitals[k], tau);
					}
					else if (atomnumber[k] == atomnum) {
						sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) * LCGTO.getSderiv(orbitals[k], orbitals[j], tau);
					}
				}
				
				H.put(j,  k, sum);
				H.put(k,  j, sum);
			}
		}
		
		DoubleMatrix G = DoubleMatrix.zeros(orbitals.length,  orbitals.length);
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {
								
								for (int l: index[a]) {
									if (l > -1) {
										for (int m: index[a]) {
											if (m > -1) {
												sum += densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
											}
										}
									}
								}
							}
						}
					}
					else {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnum]) {
									if (m > -1) {
										sum += densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
									}
								}
							}
						}
					}
				}
				
				else {
					if (atomnumber[j] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[k]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[j]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
									}
								}
							}
						}
					}
				}
				
				G.put(j,  k, sum);
				G.put(k,  j, sum);
				
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
				e += MNDODerivative.crfderivfinite(atoms[atomnum], atoms[j], tau);
			}
		}
		
		return e;

	}
	
	public static DoubleMatrix Hderiv (MNDOAtom[] atoms, DoubleMatrix densitymatrix, int atomnum, int tau) {
		int i = 0;
		
		for (MNDOAtom a: atoms) {
			i += a.getOrbitals().length;
		}
		
		MNDO6G[] orbitals = new MNDO6G[i];
		
		i = 0;
		
		int[][] index = new int[atoms.length][4];
		int[] atomnumber = new int [orbitals.length];
		int count = 0;
		int count2 = 0;
		for (MNDOAtom a: atoms) {
			count2 = 0;
			for (MNDO6G orbital: a.getOrbitals()) {
				orbitals[i] = orbital;
				index[count][count2] = i;
				atomnumber[i] = count;
				i++;
				count2++;
			}
			
			
			if (a.getZ() == 1) {
				index[count][1] = -1;
				index[count][2] = -1;
				index[count][3] = -1;
				
			}
			count++;
		}
		
		int[][] missingindex = new int[atoms.length][4 * atoms.length - 4];
		
		for (int j = 0; j < atoms.length; j++) {
			for (int k = 0; k < 4 * atoms.length - 4; k++) {
				missingindex[j][k] = -1;
			}
		}
		
		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[] {index[j][0], index[j][1], index[j][2], index[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (nums[0] != k && nums[1] != k && nums[2] != k && nums[3] != k) {
					missingindex[j][counter] = k;
					counter++;
				}
			}
		}
		
		DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {
								
								
								sum -= atoms[a].getQ() * MNDODerivative.getGderiv(orbitals[j], orbitals[k], atoms[a].s(), atoms[a].s(), tau);
							}
						}
					}
					else {
						sum -= atoms[atomnum].getQ() * MNDODerivative.getGderiv(atoms[atomnum].s(), atoms[atomnum].s(), orbitals[j], orbitals[k], tau);
					}
				}
				
				else {
					if (atomnumber[j] == atomnum) {
						sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) * LCGTO.getSderiv(orbitals[j], orbitals[k], tau);
					}
					else if (atomnumber[k] == atomnum) {
						sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) * LCGTO.getSderiv(orbitals[k], orbitals[j], tau);
					}
				}
				
				H.put(j,  k, sum);
				H.put(k,  j, sum);
			}
		}
		
		return H;
	}
	
	public static double getGderivfinite (MNDO6G a, MNDO6G b, MNDO6G c, MNDO6G d, int tau) {
		
		double orig = MNDO6G.getG(a, b, c, d);
		
		double[] newcoords = a.getCoords().clone();
		
		newcoords[tau] += 1E-8;
		
		MNDO6G anew = new MNDO6G (a, newcoords);
		MNDO6G bnew = new MNDO6G (b, newcoords);
		
		double perturbed = MNDO6G.getG(anew, bnew, c, d);
		
		return (perturbed - orig) / 1E-8;
	}
	
	public static double Hderivative (MNDOAtom[] atoms, DoubleMatrix densitymatrix, int atomnum, int tau) {
		int i = 0;
		
		for (MNDOAtom a: atoms) {
			i += a.getOrbitals().length;
		}
		
		MNDO6G[] orbitals = new MNDO6G[i];
		
		i = 0;
		
		int[][] index = new int[atoms.length][4];
		int[] atomnumber = new int [orbitals.length];
		int count = 0;
		int count2 = 0;
		for (MNDOAtom a: atoms) {
			count2 = 0;
			for (MNDO6G orbital: a.getOrbitals()) {
				orbitals[i] = orbital;
				index[count][count2] = i;
				atomnumber[i] = count;
				i++;
				count2++;
			}
			
			
			if (a.getZ() == 1) {
				index[count][1] = -1;
				index[count][2] = -1;
				index[count][3] = -1;
				
			}
			count++;
		}
		
		int[][] missingindex = new int[atoms.length][4 * atoms.length - 4];
		
		for (int j = 0; j < atoms.length; j++) {
			for (int k = 0; k < 4 * atoms.length - 4; k++) {
				missingindex[j][k] = -1;
			}
		}
		
		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[] {index[j][0], index[j][1], index[j][2], index[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (nums[0] != k && nums[1] != k && nums[2] != k && nums[3] != k) {
					missingindex[j][counter] = k;
					counter++;
				}
			}
		}
		
		DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {
								
								
								sum -= atoms[a].getQ() * MNDODerivative.getGderiv(orbitals[j], orbitals[k], atoms[a].s(), atoms[a].s(), tau);
							}
						}
					}
					else {
						sum -= atoms[atomnum].getQ() * MNDODerivative.getGderiv(atoms[atomnum].s(), atoms[atomnum].s(), orbitals[j], orbitals[k], tau);
					}
				}
				
				else {
					if (atomnumber[j] == atomnum) {
						sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) * LCGTO.getSderiv(orbitals[j], orbitals[k], tau);
					}
					else if (atomnumber[k] == atomnum) {
						sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) * LCGTO.getSderiv(orbitals[k], orbitals[j], tau);
					}
				}
				
				H.put(j,  k, sum);
				H.put(k,  j, sum);
			}
		}
		
		double e = 0;
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += densitymatrix.get(j, k) * H.get(j, k);
			}
		}
		
		return e;
	}
	
	public static DoubleMatrix Gderiv (MNDOAtom[] atoms, DoubleMatrix densitymatrix, int atomnum, int tau) {
		int i = 0;
		
		for (MNDOAtom a: atoms) {
			i += a.getOrbitals().length;
		}
		
		MNDO6G[] orbitals = new MNDO6G[i];
		
		i = 0;
		
		int[][] index = new int[atoms.length][4];
		int[] atomnumber = new int [orbitals.length];
		int count = 0;
		int count2 = 0;
		for (MNDOAtom a: atoms) {
			count2 = 0;
			for (MNDO6G orbital: a.getOrbitals()) {
				orbitals[i] = orbital;
				index[count][count2] = i;
				atomnumber[i] = count;
				i++;
				count2++;
			}
			
			
			if (a.getZ() == 1) {
				index[count][1] = -1;
				index[count][2] = -1;
				index[count][3] = -1;
				
			}
			count++;
		}
		
		int[][] missingindex = new int[atoms.length][4 * atoms.length - 4];
		
		for (int j = 0; j < atoms.length; j++) {
			for (int k = 0; k < 4 * atoms.length - 4; k++) {
				missingindex[j][k] = -1;
			}
		}
		
		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[] {index[j][0], index[j][1], index[j][2], index[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (nums[0] != k && nums[1] != k && nums[2] != k && nums[3] != k) {
					missingindex[j][counter] = k;
					counter++;
				}
			}
		}
		
		DoubleMatrix G = DoubleMatrix.zeros(orbitals.length,  orbitals.length);
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {
								
								for (int l: index[a]) {
									if (l > -1) {
										for (int m: index[a]) {
											if (m > -1) {
												sum += densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
											}
										}
									}
								}
							}
						}
					}
					else {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnum]) {
									if (m > -1) {
										sum += densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
									}
								}
							}
						}
					}
				}
				
				else {
					if (atomnumber[j] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[k]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[j]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
									}
								}
							}
						}
					}
				}
				
				G.put(j,  k, sum);
				G.put(k,  j, sum);
				
			}
		}
		
		return G;
	}
	
	public static double Gderivative (MNDOAtom[] atoms, DoubleMatrix densitymatrix, int atomnum, int tau) {
		int i = 0;
		
		for (MNDOAtom a: atoms) {
			i += a.getOrbitals().length;
		}
		
		MNDO6G[] orbitals = new MNDO6G[i];
		
		i = 0;
		
		int[][] index = new int[atoms.length][4];
		int[] atomnumber = new int [orbitals.length];
		int count = 0;
		int count2 = 0;
		for (MNDOAtom a: atoms) {
			count2 = 0;
			for (MNDO6G orbital: a.getOrbitals()) {
				orbitals[i] = orbital;
				index[count][count2] = i;
				atomnumber[i] = count;
				i++;
				count2++;
			}
			
			
			if (a.getZ() == 1) {
				index[count][1] = -1;
				index[count][2] = -1;
				index[count][3] = -1;
				
			}
			count++;
		}
		
		int[][] missingindex = new int[atoms.length][4 * atoms.length - 4];
		
		for (int j = 0; j < atoms.length; j++) {
			for (int k = 0; k < 4 * atoms.length - 4; k++) {
				missingindex[j][k] = -1;
			}
		}
		
		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[] {index[j][0], index[j][1], index[j][2], index[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (nums[0] != k && nums[1] != k && nums[2] != k && nums[3] != k) {
					missingindex[j][counter] = k;
					counter++;
				}
			}
		}
		
		DoubleMatrix G = DoubleMatrix.zeros(orbitals.length,  orbitals.length);
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {
								
								for (int l: index[a]) {
									if (l > -1) {
										for (int m: index[a]) {
											if (m > -1) {
												sum += densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
											}
										}
									}
								}
							}
						}
					}
					else {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnum]) {
									if (m > -1) {
										sum += densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
									}
								}
							}
						}
					}
				}
				
				else {
					if (atomnumber[j] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[k]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[j]]) {
									if (m > -1) {
										sum -= 0.5 * densitymatrix.get(l, m) * MNDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
									}
								}
							}
						}
					}
				}
				
				G.put(j,  k, sum);
				G.put(k,  j, sum);
				
			}
		}
		
		double e = 0;
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += densitymatrix.get(j, k) * G.get(j, k);
			}
		}
		
		return e;
	}
	
	public static double crfderivative (MNDOAtom[] atoms, int atomnum, int tau) {
		
		double e = 0;
		for (int j = 0; j < atoms.length; j++) {
			if (j != atomnum) {
				e += MNDODerivative.crfderivfinite(atoms[atomnum], atoms[j], tau);
			}
		}
		
		return e;
		
	}
	
	
	
	public static double gradientfinite (MNDOAtom[] atoms, int charge, int atomnum, int tau) {
		MNDOSolution s = new MNDOSolution (atoms, charge);
		
		MNDOAtom[] perturbed = new MNDOAtom [atoms.length];
		
		for (int i = 0; i < atoms.length; i++) {
			perturbed[i] = new MNDOAtom(atoms[i]);
			
			if (i == atomnum) {
				perturbed[i].getCoordinates()[tau] = perturbed[i].getCoordinates()[tau] + 1E-8;
			}
		}
		
		MNDOSolution sprime = new MNDOSolution (perturbed, charge);
		
		return 1E8 * (sprime.energy-s.energy);
	}
	
	public static double gradientunrestricted (MNDOAtom[] atoms, DoubleMatrix alphadensity, DoubleMatrix betadensity, int atomnum, int tau) {
		int i = 0;
		
		for (MNDOAtom a: atoms) {
			i += a.getOrbitals().length;
		}
		
		MNDO6G[] orbitals = new MNDO6G[i];
		
		i = 0;
		
		int[][] index = new int[atoms.length][4];
		int[] atomnumber = new int [orbitals.length];
		int count = 0;
		int count2 = 0;
		for (MNDOAtom a: atoms) {
			count2 = 0;
			for (MNDO6G orbital: a.getOrbitals()) {
				orbitals[i] = orbital;
				index[count][count2] = i;
				atomnumber[i] = count;
				i++;
				count2++;
			}
			
			
			if (a.getZ() == 1) {
				index[count][1] = -1;
				index[count][2] = -1;
				index[count][3] = -1;
				
			}
			count++;
		}
		
		int[][] missingindex = new int[atoms.length][4 * atoms.length - 4];
		
		for (int j = 0; j < atoms.length; j++) {
			for (int k = 0; k < 4 * atoms.length - 4; k++) {
				missingindex[j][k] = -1;
			}
		}
		
		for (int j = 0; j < atoms.length; j++) {
			int[] nums = new int[] {index[j][0], index[j][1], index[j][2], index[j][3]};
			int counter = 0;
			for (int k = 0; k < orbitals.length; k++) {
				if (nums[0] != k && nums[1] != k && nums[2] != k && nums[3] != k) {
					missingindex[j][counter] = k;
					counter++;
				}
			}
		}
		
		DoubleMatrix H = DoubleMatrix.zeros(orbitals.length, orbitals.length);
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				if (atomnumber[j] == atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int a = 0; a < atoms.length; a++) {
							if (a != atomnum) {
								
								
								sum -= atoms[a].getQ() * MNDODerivative.getGderiv(orbitals[j], orbitals[k], atoms[a].s(), atoms[a].s(), tau);
							}
						}
					}
					else {
						sum -= atoms[atomnum].getQ() * MNDODerivative.getGderiv(atoms[atomnum].s(), atoms[atomnum].s(), orbitals[j], orbitals[k], tau);
					}
				}
				
				else {
					if (atomnumber[j] == atomnum) {
						sum += 0.5 * (orbitals[j].beta() + orbitals[k].beta()) * LCGTO.getSderiv(orbitals[j], orbitals[k], tau);
					}
					else if (atomnumber[k] == atomnum) {
						sum += 0.5 * (orbitals[k].beta() + orbitals[j].beta()) * LCGTO.getSderiv(orbitals[k], orbitals[j], tau);
					}
				}
				
				H.put(j,  k, sum);
				H.put(k,  j, sum);
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
								
								for (int l: index[a]) {
									if (l > -1) {
										for (int m: index[a]) {
											if (m > -1) {
												sum += (alphadensity.get(l, m) + betadensity.get(l, m)) * MNDODerivative.getGderiv(orbitals[j], orbitals[k], orbitals[l], orbitals[m], tau);
											}
										}
									}
								}
							}
						}
					}
					else {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnum]) {
									if (m > -1) {
										sum += (alphadensity.get(l, m) + betadensity.get(l, m)) * MNDODerivative.getGderiv(orbitals[l], orbitals[m], orbitals[j], orbitals[k], tau);
									}
								}
							}
						}
					}
				}
				
				J.put(j,  k, sum);
				J.put(k,  j, sum);
				
			}
		}
		
		DoubleMatrix Ka = DoubleMatrix.zeros(orbitals.length, orbitals.length);
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				
				if (atomnumber[j] != atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[k]]) {
									if (m > -1) {
										sum -= alphadensity.get(l, m) * MNDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[j]]) {
									if (m > -1) {
										sum -= alphadensity.get(l, m) * MNDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
									}
								}
							}
						}
					}
				}
				
				Ka.put(j,  k, sum);
				Ka.put(k,  j, sum);
				
			}
		}
		
		DoubleMatrix Kb = DoubleMatrix.zeros(orbitals.length, orbitals.length);
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				double sum = 0;
				
				if (atomnumber[j] != atomnumber[k]) {
					if (atomnumber[j] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[k]]) {
									if (m > -1) {
										sum -= betadensity.get(l, m) * MNDODerivative.getGderiv(orbitals[j], orbitals[l], orbitals[k], orbitals[m], tau);
									}
								}
							}
						}
					}
					else if (atomnumber[k] == atomnum) {
						for (int l: index[atomnum]) {
							if (l > -1) {
								for (int m: index[atomnumber[j]]) {
									if (m > -1) {
										sum -= betadensity.get(l, m) * MNDODerivative.getGderiv(orbitals[k], orbitals[l], orbitals[j], orbitals[m], tau);
									}
								}
							}
						}
					}
				}
				
				Kb.put(j,  k, sum);
				Kb.put(k,  j, sum);
				
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
				e += MNDODerivative.crfderivfinite(atoms[atomnum], atoms[j], tau);
			}
		}
		
		return e;
	}

}
