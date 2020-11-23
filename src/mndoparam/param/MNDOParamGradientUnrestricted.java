package mndoparam.param;

import java.util.Arrays;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDOGeometryOptimizationUnrestricted;
import mndoparam.mndo.MNDOSolutionUnrestricted;

public class MNDOParamGradientUnrestricted {
	
	private static final double lambda = 1E-7;
	
	
	public MNDOSolutionUnrestricted s;
	
	private MNDOSolutionUnrestricted sprime;
	
	public MNDOParamErrorFunctionUnrestricted e, eprime;
	
	private MNDOAtom[] atoms, perturbed;
	
	public String str;
	
	private int Z, paramnum;
	
	public MNDOParamGradientUnrestricted (MNDOAtom[] atoms, int charge, int mult, int Z, int paramnum) {
		
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
		
		s = new MNDOGeometryOptimizationUnrestricted(atoms, charge, mult).s;
		
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
		
		sprime = new MNDOGeometryOptimizationUnrestricted(perturbed, charge, mult).s;
		
		
	}
	
	public MNDOParamGradientUnrestricted (MNDOAtom[] atoms, int charge, int mult, int Z, int paramnum, MNDOGeometryOptimizationUnrestricted s) {
		
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
		
		sprime = new MNDOSolutionUnrestricted (perturbed, charge, mult);
		
		//sprime = new MNDOParamGeometryOptimizationUnrestricted (perturbed, charge, mult).s;
	}
	
	public MNDOParamGradientUnrestricted (MNDOAtom[] atoms, int charge, int mult, int Z, int paramnum, MNDOSolutionUnrestricted s) {
		
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
		
		sprime = new MNDOSolutionUnrestricted (perturbed, charge, mult);
		
		//sprime = new MNDOParamGeometryOptimizationUnrestricted (perturbed, charge, mult).s;
	}
	
	public void constructErrors (double refHeat) {
		e = new MNDOParamErrorFunctionUnrestricted (atoms, s, refHeat);
		eprime = new MNDOParamErrorFunctionUnrestricted (perturbed, sprime, refHeat);
	}
	
	public void addDipoleError (double ref) {
		
		e.AddDipoleError(ref);
		eprime.AddDipoleError(ref);
	}
	
	public void addIEError (double ref) {
		e.AddIEError(ref);
		eprime.AddIEError(ref);
	}
	
	public void createExpGeom (MNDOAtom[] expatoms, MNDOSolutionUnrestricted expsoln) {
		e.createExpGeom(expatoms, expsoln);
		
		MNDOAtom[] perturbed = new MNDOAtom [expatoms.length];
		
		for (int i = 0; i < expatoms.length; i++) {
			
			perturbed[i] = new MNDOAtom (expatoms[i]);
			
			if (expatoms[i].getZ() == Z) {
				double[] params = expatoms[i].getParams();
				
				//System.err.println(Arrays.toString(params));
				
				params[paramnum] = params[paramnum] + lambda;
				
				//System.err.println(Arrays.toString(params));
				
				perturbed[i] = new MNDOAtom (expatoms[i].getCoordinates(), Z, params);
			}
			
		}
		eprime.createExpGeom(perturbed, new MNDOSolutionUnrestricted (perturbed, expsoln.charge, expsoln.multiplicity));
		
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
