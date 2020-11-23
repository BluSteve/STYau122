package mndoparam.mndo;

import java.util.Arrays;

import org.jblas.DoubleMatrix;
import org.jblas.Solve;

public class MNDOGeometryOptimization {
	private MNDOAtom[] atoms;
	
	private double refEnergy;
	
	private int counter;
	
	public double IE, dipole;

	public double heat;
	
	private DoubleMatrix densitymatrix;
	
	public MNDOSolution s;
	
	public MNDOGeometryOptimization (MNDOAtom[] atoms, int charge) {
		
		this.atoms = atoms;
		
		//System.out.println ("Initial coordinates:\n");
		
		for (MNDOAtom a: atoms) {
			//System.out.println (Arrays.toString(a.getCoordinates()));
		}
		
		
		//System.out.println ("\nObtaining reference energy...");
		
		s = new MNDOSolution (atoms, charge);
		
		refEnergy = s.energy;
		
		densitymatrix = s.densitymatrix();
		
		System.out.println ("\nCurrent heat of formation: " + s.hf + "kcal/mol");
		System.out.println ("Current HOMO energy: " + s.homo + " eV");
		//System.out.println ("Current LUMO energy: " + s.lumo + " eV");
		
		System.out.println ("-----------------------------------------------");
		
		DoubleMatrix gradient = new DoubleMatrix (atoms.length * 3, 1);
		
		int index = 0;
		
		for (int a = 0; a < atoms.length; a++) {
			for (int i =0; i < 3; i++) {
				gradient.put(index,  0, derivative (a, i));
				index++;
			}
		}
		
		//System.out.println ("gradient: " + gradient);
		
		double sum = 0;
		
		for (int i = 0; i < gradient.length; i++) {
			sum += gradient.get(i) * gradient.get(i);
		}
		
		sum = Math.sqrt(sum);
		
		
		DoubleMatrix B = DoubleMatrix.eye(atoms.length * 3);
		
		DoubleMatrix searchdir = new DoubleMatrix (atoms.length * 3, 1);
		
		
		searchdir = gradient.mul(-1/sum);
		
		DoubleMatrix oldgrad = DoubleMatrix.zeros(atoms.length * 3, 1);
		
		double energy = refEnergy - 1;
		

		double scale = 0.1;
		
		int count = 0;
	

		
		
		
		
		counter = 0;
		
		while (mag (gradient) > 0.04) {
			
			System.out.println ("Gradient: " + mag(gradient));
			
			count = 0;
			
			scale = 0.1;
			
			double summ = 0;
			
			refEnergy = 0;
			
			
			
			while (Math.abs(energy - refEnergy)> 1E-9) {
				
				
				refEnergy = energy;
				
				count = 0;
				
				for (MNDOAtom a: atoms) {
					for (int i = 0; i < 3; i++) {
						a.getCoordinates()[i] = Math.round((a.getCoordinates()[i] + scale * searchdir.get(count)) * 1000000000)/1000000000.0;
						count++;
					}
				}
				
				summ += scale;
				//System.out.println ("New coordinates:");
				
				for (MNDOAtom a: atoms) {
					//System.out.println (Arrays.toString(a.getCoordinates()));
				}
				
				//System.out.println ("\nObtaining reference energy...");
				
				s = new MNDOSolution (atoms, charge);
				
				//System.out.println ("Energy: " + s.energy + "eV\n");
				System.out.println ("\nCurrent heat of formation: " + s.hf + "kcal/mol");
				System.out.println ("Current HOMO energy: " + s.homo + " eV");
				//System.out.println ("Current LUMO energy: " + s.lumo + " eV");
				
				System.out.println ("-----------------------------------------------");
				
				energy = s.energy;
				
				
				
				if (energy > refEnergy) {
					scale *= -0.5;
				}
				
			}
			
			densitymatrix = s.densitymatrix();
			
			refEnergy = energy;
			
			oldgrad = gradient.dup();
			
			gradient = new DoubleMatrix (atoms.length * 3, 1);
			
			index = 0;
			
			for (int a = 0; a < atoms.length; a++) {
				for (int i =0; i < 3; i++) {
					gradient.put(index,  0, derivative (a, i));
					index++;
				}
			}
			
			
			sum = 0;
			
			DoubleMatrix y = gradient.sub(oldgrad);
			
			try {
				B = getb (B, y, searchdir);
				
				
				searchdir = Solve.pinv(B).mmul(gradient).mmul(-1);
			}
			catch (Exception e) {
				B = DoubleMatrix.eye(atoms.length * 3);
				
				searchdir = Solve.pinv(B).mmul(gradient).mmul(-1);
				
			}
			
			for (int i = 0; i < gradient.length; i++) {
				sum += searchdir.get(i) * searchdir.get(i);
			}
			
			sum = Math.sqrt(sum);
			
			searchdir = searchdir.mmul(0.1/sum);
		
			
			
			counter++;
		}
		
		System.out.println ("FINAL:");
		
		for (MNDOAtom a: atoms) {
			//System.out.println (Arrays.toString(a.getCoordinates()));
		}
		//System.out.println();
		
		s = new MNDOSolution (atoms, charge);
		
		//System.out.println ("Eigenvalues: " + s.E + "\n");
		
		System.out.println ("\nHeat of formation: " + s.hf + "kcal/mol");
		System.out.println ("HOMO energy: " + s.homo + " eV");
		//System.out.println ("LUMO energy: " + s.lumo + " eV");
		//System.out.println ("Total energy: " + s.energy + "eV");
		
		//System.out.println();
		
		//System.out.println ("===DIPOLE MOMENT EVALUATION===");
		
		//System.out.println ("point charge dipole contribution: " + Arrays.toString(s.chargedip));
		
		//System.out.println ("hybrid dipole contribution: " + Arrays.toString(s.hybridip));
		
		//System.out.println ("sum: " + Arrays.toString(s.dipoletot));
		
		//System.out.println ("dipole moment: " + s.dipole + " debyes");
		
		//System.out.println ("===END OF DIPOLE MOMENT DATA===");
		
		//System.out.println ("\n" + counter + " iterations");
		
		this.dipole = s.dipole;
		this.heat = s.hf;
		this.IE = -s.homo;
	}
	
	private DoubleMatrix getb (DoubleMatrix B, DoubleMatrix y, DoubleMatrix searchdir) {
		
		double a = 1 / y.transpose().mmul(searchdir).get(0);
		
		double b = searchdir.transpose().mmul(B).mmul(searchdir).get(0);
		
		DoubleMatrix m2 = B.mmul(searchdir).mmul(searchdir.transpose()).mmul(B.transpose()).mmul(b);
		
		DoubleMatrix m1 = y.mmul(y.transpose()).mmul(a);
		
		
		return B.add(m1).sub(m2);
	}
	
	private double derivative (int i, int j) {
		//System.out.println ("Evaluating gradient analytically; give us a moment...");
		return MNDODerivative.gradient(atoms, densitymatrix, i, j);
	}
	
	private double mag (DoubleMatrix gradient) {
		
		double sum = 0;
		for (int i = 0; i < gradient.length; i++) {
			sum += gradient.get(i) * gradient.get(i);
		}
		
		return Math.sqrt(sum);
	}
	
	public MNDOAtom[] getAtoms() {
		return this.atoms;
	}
}
