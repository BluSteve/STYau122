package mndoparam.param;

import java.util.Arrays;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDOGeometryOptimization;
import mndoparam.mndo.MNDOSolution;

public class MNDOParamGradient {
	
	public static final double lambda = 1E-7;
	
	public MNDOSolution s;
	
	private MNDOSolution sprime;
	
	public MNDOParamErrorFunction e, eprime;
	
	private MNDOAtom[] atoms;

	public MNDOAtom[] perturbed;
	
	private int Z, paramnum;
	
	public String str;
	
	public MNDOParamGradient (MNDOAtom[] atoms, int charge, int Z, int paramnum) {
		
		this.Z = Z;
		this.paramnum = paramnum;
		str = "";
		
		switch (Z) {
		case 1:
			System.err.print ("H ");
			str += "H ";
			break;
		case 6:
			System.err.print ("C ");
			str += "C ";
			break;
		case 7:
			System.err.print ("N ");
			str += "N ";
			break;
		}
		
		switch (paramnum) {
		case 0:
			System.err.println ("ALPHA");
			str += "ALPHA: ";
			break;
		case 1:
			System.err.println ("BETAS");
			str += "BETAS: ";
			break;
		case 2:
			System.err.println ("BETAP");
			str += "BETAP: ";
			break;
		case 3:
			System.err.println ("USS");
			str += "USS: ";
			break;
		case 4:
			System.err.println ("UPP");
			str += "UPP: ";
			break;
		case 5:
			System.err.println ("ZETAS");
			str += "ZETAS: ";
			break;
		case 6:
			System.err.println ("ZETAP");
			str += "ZETAP: ";
			break;
		case 7:
			System.err.println ("EISOL");
			str += "EISOL: ";
			break;
		default:
			System.err.println ("DO NOT TOUCH");
			System.exit(0);
			
		}
		
		this.atoms = atoms;
		
		s = new MNDOGeometryOptimization(atoms, charge).s;
		
		perturbed = new MNDOAtom [atoms.length];
		
		for (int i = 0; i < atoms.length; i++) {
			
			perturbed[i] = new MNDOAtom (atoms[i]);
			
			if (atoms[i].getZ() == Z) {
				double[] params = atoms[i].getParams();
				
				//System.err.println(Arrays.toString(params));
				
				params[paramnum] = params[paramnum] + lambda;
				
				//System.err.println(Arrays.toString(params));
				
				perturbed[i] = new MNDOAtom (atoms[i].getCoordinates(), Z, params);
			}
		}
		
		sprime = new MNDOGeometryOptimization(perturbed, charge).s;
		
		
	}
	
	public MNDOParamGradient (MNDOAtom[] atoms, int charge, int Z, int paramnum, MNDOGeometryOptimization s) {
		this.Z = Z;
		this.paramnum = paramnum;
		str = "";
		
		switch (Z) {
		case 1:
			System.err.print ("H ");
			str += "H ";
			break;
		case 6:
			System.err.print ("C ");
			str += "C ";
			break;
		case 7:
			System.err.print ("N ");
			str += "N ";
			break;
		}
		
		switch (paramnum) {
		case 0:
			System.err.println ("ALPHA");
			str += "ALPHA: ";
			break;
		case 1:
			System.err.println ("BETAS");
			str += "BETAS: ";
			break;
		case 2:
			System.err.println ("BETAP");
			str += "BETAP: ";
			break;
		case 3:
			System.err.println ("USS");
			str += "USS: ";
			break;
		case 4:
			System.err.println ("UPP");
			str += "UPP: ";
			break;
		case 5:
			System.err.println ("ZETAS");
			str += "ZETAS: ";
			break;
		case 6:
			System.err.println ("ZETAP");
			str += "ZETAP: ";
			break;
		case 7:
			System.err.println ("EISOL");
			str += "EISOL: ";
			break;
		default:
			System.err.println ("DO NOT TOUCH");
			System.exit(0);
			
		}
		
		this.atoms = atoms;
		
		this.s = s.s;
		
		perturbed = new MNDOAtom [atoms.length];
		
		for (int i = 0; i < atoms.length; i++) {
			
			perturbed[i] = new MNDOAtom (atoms[i]);
			
			if (atoms[i].getZ() == Z) {
				double[] params = atoms[i].getParams();
				
				//System.err.println(Arrays.toString(params));
				
				params[paramnum] = params[paramnum] + lambda;
				
				//System.err.println(Arrays.toString(params));
				
				perturbed[i] = new MNDOAtom (atoms[i].getCoordinates(), Z, params);
			}
		}
		
		sprime = new MNDOSolution(perturbed, charge);
		
		//sprime = new MNDOParamGeometryOptimization (perturbed, charge).s;
	}
	
	public MNDOParamGradient (MNDOAtom[] atoms, int charge, int Z, int paramnum, MNDOSolution s) {
		this.Z = Z;
		this.paramnum = paramnum;
		str = "";
		
		switch (Z) {
		case 1:
			System.err.print ("H ");
			str += "H ";
			break;
		case 6:
			System.err.print ("C ");
			str += "C ";
			break;
		case 7:
			System.err.print ("N ");
			str += "N ";
			break;
		}
		
		switch (paramnum) {
		case 0:
			System.err.println ("ALPHA");
			str += "ALPHA: ";
			break;
		case 1:
			System.err.println ("BETAS");
			str += "BETAS: ";
			break;
		case 2:
			System.err.println ("BETAP");
			str += "BETAP: ";
			break;
		case 3:
			System.err.println ("USS");
			str += "USS: ";
			break;
		case 4:
			System.err.println ("UPP");
			str += "UPP: ";
			break;
		case 5:
			System.err.println ("ZETAS");
			str += "ZETAS: ";
			break;
		case 6:
			System.err.println ("ZETAP");
			str += "ZETAP: ";
			break;
		case 7:
			System.err.println ("EISOL");
			str += "EISOL: ";
			break;
		default:
			System.err.println ("DO NOT TOUCH");
			System.exit(0);
			
		}
		
		this.atoms = atoms;
		
		this.s = s;
		
		perturbed = new MNDOAtom [atoms.length];
		
		for (int i = 0; i < atoms.length; i++) {
			
			perturbed[i] = new MNDOAtom (atoms[i]);
			
			if (atoms[i].getZ() == Z) {
				double[] params = atoms[i].getParams();
				
				//System.err.println(Arrays.toString(params));
				
				params[paramnum] = params[paramnum] + lambda;
				
				//System.err.println(Arrays.toString(params));
				
				perturbed[i] = new MNDOAtom (atoms[i].getCoordinates(), Z, params);
			}
		}
		
		sprime = new MNDOSolution(perturbed, charge);
		
		//sprime = new MNDOParamGeometryOptimization (perturbed, charge).s;
	}
	
	public void constructErrors (double refHeat) {
		e = new MNDOParamErrorFunction (atoms, s, refHeat);
		eprime = new MNDOParamErrorFunction (perturbed, sprime, refHeat);
	}
	
	public void addDipoleError (double ref) {
		
		e.AddDipoleError(ref);
		eprime.AddDipoleError(ref);
	}
	
	public void addIEError (double ref) {
		e.AddIEError(ref);
		eprime.AddIEError(ref);
	}
	
	public void createExpGeom (MNDOAtom[] expatoms, MNDOSolution expsoln) {
		e.createExpGeom(expatoms, expsoln);
		
		MNDOAtom[] perturbed = new MNDOAtom [expatoms.length];
		
		for (int i = 0; i < expatoms.length; i++) {
			
			perturbed[i] = new MNDOAtom (expatoms[i]);
			
			//System.err.println (Z);
			
			if (expatoms[i].getZ() == Z) {
				double[] params = expatoms[i].getParams();
				
				//System.err.println(Arrays.toString(params));
				
				params[paramnum] = params[paramnum] + lambda;
				
				//System.err.println(Arrays.toString(params));
				
				perturbed[i] = new MNDOAtom (expatoms[i], params);
			}
		}
		eprime.createExpGeom(perturbed, new MNDOSolution (perturbed, expsoln.charge));
	}
	
	public void addGeomError () {
		e.AddGeomError();
		eprime.AddGeomError();
	}
	
	public void addBondError (int atom1, int atom2, double ref) {
		e.AddBondError(atom1, atom2, ref);
		eprime.AddBondError(atom1, atom2, ref);
	}
	
	public void addAngleError (int atom1, int atom2, int atom3, double ref) {
		e.AddAngleError(atom1, atom2, atom3, ref);
		eprime.AddAngleError(atom1, atom2, atom3, ref);
	}
	
	public double gradient() {
		return (eprime.constructErrorFunction() - e.constructErrorFunction()) / lambda;
	}

}
