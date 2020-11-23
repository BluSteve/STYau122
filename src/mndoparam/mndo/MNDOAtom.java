package mndoparam.mndo;

import scf.Atom;
import scf.GTO;
import scf.LCGTO;

public class MNDOAtom extends Atom{
	
	protected double alpha, betas, betap, Uss, Upp, zetas, zetap, eisol, gss, gsp, hsp, gpp, gp2, p0, p1, p2, D1, D2;
	
	private MNDO6G[] orbitals;

	public MNDOAtom(double[] coordinates, int Z, double alpha, double betas, double betap, double Uss, double Upp, double zetas, double zetap, double eisol, double gss, double gsp, double hsp, double gpp, double gp2) {
		super(new MNDO6G[] {}, coordinates, Z);
		
		this.alpha = alpha;
		this.betas = betas;
		this.betap = betap;
		this.Uss = Uss;
		this.Upp = Upp;
		this.zetas = zetas;
		this.zetap = zetap;
		this.eisol = eisol;
		this.gss = gss;
		this.gsp = gsp;
		this.hsp = hsp;
		this.gpp = gpp;
		this.gp2 = gp2;
		this.p0 = p0();
		this.D1 = D1();
		this.D2 = D2();
		this.p1 = p1();
		this.p2 = p2();
		
		//System.out.println ("Derived parameters: " + p0 + ", " + p1 + ", " + p2 + ", " + D1 + ", " + D2);
		
		this.orbitals = getOrbitals (coordinates, Z);
		
		super.orbitals = getOrbitals (coordinates, Z);
	}
	
	public MNDO6G[] getOrbitals (double[] coordinates, int Z) {//initialises OM2-3G basis functions
		switch (Z) {
		case 1:
			return new MNDO6G[] {new MNDO6G(0, 0, 0, coordinates, 1, Z, this, zetas, betas, Uss, p0, 0, 0, 0, 0, gss, 0, 0, 0, 0)};
		case 6:
		case 7:
		case 8:
		case 9:
			return new MNDO6G[] {new MNDO6G(0, 0, 0, coordinates, 2, Z, this, zetas, betas, Uss, p0, p1, p2, D1, D2, gss, gsp, hsp, gpp, gp2), 
								new MNDO6G(1, 0, 0, coordinates, 2, Z, this, zetap, betap, Upp, p0, p1, p2, D1, D2, gss, gsp, hsp, gpp, gp2),
								new MNDO6G(0, 1, 0, coordinates, 2, Z, this, zetap, betap, Upp, p0, p1, p2, D1, D2, gss, gsp, hsp, gpp, gp2),
								new MNDO6G(0, 0, 1, coordinates, 2, Z, this, zetap, betap, Upp, p0, p1, p2, D1, D2, gss, gsp, hsp, gpp, gp2)};
		default:
			return new MNDO6G[] {};
		
		}
	}
	
	public double mass () {
		switch (getZ()) {
		
		case 1:
			return 1.00797;
		case 6:
			return 12.011;
		case 7:
			return 14.0067;
		case 8:
			return 15.9994;
		case 9:
			return 18.998403;
		}
		
		return 0;
	}
	
	public double[] getParams() {
		return new double[] {alpha, betas, betap, Uss, Upp, zetas, zetap, eisol, gss, gsp, hsp, gpp, gp2}.clone();
	}
	
	public MNDOAtom (double[] coords, MNDOAtom a) {
		this(coords.clone(), a.getZ(), a.alpha, a.betas, a.betap, a.Uss, a.Upp, a.zetas, a.zetap, a.eisol, a.gss, a.gsp, a.hsp, a.gpp, a.gp2);
	}
	
	public MNDOAtom (double[] coords, int Z, double[] params) {
		this (coords.clone(), Z, params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], params[10], params[11], params[12]);
	}
	
	public MNDOAtom (MNDOAtom a) {
		this(a.getCoordinates().clone(), a.getZ(), a.alpha, a.betas, a.betap, a.Uss, a.Upp, a.zetas, a.zetap, a.eisol, a.gss, a.gsp, a.hsp, a.gpp, a.gp2);
	}
	
	public MNDOAtom (MNDOAtom a, double[] params) {
		this(a.getCoordinates().clone(), a.getZ(), params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9], params[10], params[11], params[12]);
	}
	
	public MNDO6G[] getOrbitals() {
		return this.orbitals;
	}
	
	public MNDO6G s() {
		return this.orbitals[0];
	}
	
	public double V (MNDO6G a, MNDO6G b) {
		
		return -this.getQ() * MNDO6G.getG(a, b, this.s(), this.s());
	}
	
	public static double crf (MNDOAtom a, MNDOAtom b) {
		double f;
		double R = GTO.R(a.getCoordinates(), b.getCoordinates()) / 1.88973;
		if ((a.getZ() == 7 || a.getZ() == 8) && b.getZ() == 1) {
			f = 1 + R * Math.exp(-a.alpha * R) + Math.exp(-b.alpha * R);
		}
		else if ((b.getZ() == 7 || b.getZ() == 8) && a.getZ() == 1) {
			f = 1 + R * Math.exp(-b.alpha * R) + Math.exp(-a.alpha * R);
		}
		else {
			f = 1 + Math.exp(-b.alpha * R) + Math.exp(-a.alpha * R);
		}
		
		return f * a.getQ() * b.getQ() * MNDO6G.getG(a.s(), a.s(), b.s(), b.s());
	}
	
	public static double crfderiv (MNDOAtom a, MNDOAtom b, int tau) {
		double f;
		double fprime;
		double R = GTO.R(a.getCoordinates(), b.getCoordinates());
		if ((a.getZ() == 7 || a.getZ() == 8) && b.getZ() == 1) {
			f = 1 + R / 1.88973 * Math.exp(-a.alpha * R / 1.88973) + Math.exp(-b.alpha * R / 1.88973);
			
			fprime = (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * 1.88973) * Math.exp(-a.alpha * R / 1.88973)
				   - a.alpha * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (1.88973 * 1.88973) * Math.exp(-a.alpha * R / 1.88973) 
				   - b.alpha / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-b.alpha * R / 1.88973);
		}
		else if ((b.getZ() == 7 || b.getZ() == 8) && a.getZ() == 1) {
			f = 1 + R / 1.88973 * Math.exp(-b.alpha * R / 1.88973) + Math.exp(-a.alpha * R / 1.88973);
			
			fprime = (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * 1.88973) * Math.exp(-b.alpha * R / 1.88973)
				   - b.alpha * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (1.88973 * 1.88973) * Math.exp(-b.alpha * R / 1.88973) 
				   - a.alpha / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-a.alpha * R / 1.88973);
		}
		else {
			f = 1 + Math.exp(-b.alpha * R / 1.88973) + Math.exp(-a.alpha * R / 1.88973);
			
			fprime = - b.alpha / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-b.alpha * R / 1.88973)
				     - a.alpha / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-a.alpha * R / 1.88973);
		}
		
		return fprime * a.getQ() * b.getQ() * MNDO6G.getG(a.s(), a.s(), b.s(), b.s()) + f * a.getQ() * b.getQ() * MNDODerivative.getGderiv(a.s(), a.s(), b.s(), b.s(), tau);
	}
	
	public double heat () {
		double conversion = 4.3363E-2;
		switch (getZ()) {
		case 1:
			return 52.10200000 * conversion;
		case 6:
			return 170.89000000 * conversion;
		case 7:
			return 113.00000000 * conversion;
		case 8:
			return 59.55900000 * conversion;
		case 9:
			return 18.89000000 * conversion; 
		}
		return 0;
	}
	
	public double eisol() {
		return eisol;
	}
	
	private double p0() {
		return 27.2114/ (2 * gss);
	}
	
	private double D1() {
		return (2 * shell + 1) / Math.sqrt(3) * Math.pow(4 * zetas * zetap, shell + 0.5) / Math.pow(zetas + zetap, 2 * shell + 2);
	}
	
	private double D2() {
		return 1/ zetap * Math.sqrt((2 * shell + 1) * (2 * shell + 2) / 20.0);
	}
	
	private double p1() {
		double guess = 0;
		
		double newguess = 0.5 * Math.pow(D1 * D1 * 27.2114 / (hsp), 1.0/3);
		
		while (Math.abs(guess - newguess) > 1E-12) {
			
			guess = newguess;
			double f = 1/guess - 1/Math.sqrt(guess * guess + D1 * D1) - 4 * hsp / 27.2114;
			double fprime = -1/(guess * guess) + guess / Math.pow(guess * guess + D1 * D1, 1.5);
			
			newguess = guess - f/fprime;
		}
		return newguess;
	}
	
	private double p2() {
		double guess = 0;
		double newguess = 0.5;
		
		while (Math.abs(guess - newguess) > 1E-12) {
			guess = newguess;
			double f = 1/guess + 1/Math.sqrt(guess * guess + 2 * D2 * D2) - 2/Math.sqrt(guess * guess + D2 * D2) - 4 * (gpp - gp2) / 27.2114;
			double fprime = -1/(guess * guess) - guess / Math.pow(guess * guess + 2 * D2 * D2, 1.5) + 2 * guess / Math.pow(guess * guess + D2 * D2, 1.5);
			
			newguess = guess - f/fprime;
		}
		return newguess;
	}
	
	

}
