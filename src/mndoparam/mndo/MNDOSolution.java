package mndoparam.mndo;

import java.util.ArrayList;
import java.util.Arrays;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import runcycle.MoleculeRun;


public class MNDOSolution {
	
	private MNDO6G[] orbitals;
	
	private MNDOAtom[] atoms;
	
	public double energy;
	
	private DoubleMatrix densitymatrix;
	
	public DoubleMatrix H, C, F, G, E;//H - core matrix, G = 2-electron matrix, F = fock matrix, C = coeffecient matrix (transposed for easier reading), E = eigenvalues 
	
	private int NElectrons;
	
	public double hf, homo, lumo;
	
	public double[] chargedip, hybridip, dipoletot;
	
	public double dipole;
	
	public int charge;
	
	private ArrayList<Double> integralarray;//my lazy implementation of integral storage
	
	public MNDOSolution (MNDOAtom[] atoms, int charge) {
		
		this.charge = charge;
		
		double damp = 0.8;
		
		int n = 0;
		
		for (MNDOAtom a: atoms) {
			n += a.getQ();
		}
		
		n -= charge;
		
		this.NElectrons = n;
		this.atoms = atoms;
		
		int i = 0;
		
		for (MNDOAtom a: atoms) {
			i += a.getOrbitals().length;
		}
		
		orbitals = new MNDO6G[i];
		
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
		
		H = new DoubleMatrix (orbitals.length, orbitals.length);
		
		//filling up the core matrix in accordance with MNDO formalism
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (k == j) {
					double Huu = orbitals[j].U();
					for (MNDOAtom a: atoms) {
						if (!Arrays.equals(a.getCoordinates(), orbitals[j].getCoords())) {
							Huu += a.V(orbitals[j], orbitals[k]); 
						}
					}
					H.put(j, k, Huu);
				}
				
				else if (Arrays.equals(orbitals[j].getCoords(), orbitals[k].getCoords())) {
					double Huv = 0;
					for (MNDOAtom a: atoms) {
						if (!Arrays.equals(a.getCoordinates(), orbitals[j].getCoords())) {
							Huv += a.V(orbitals[j], orbitals[k]); 
						}
					}
					H.put(j, k, Huv);
					H.put(k, j, Huv);
				}
				
				else {
					double Huk = MNDO6G.beta(orbitals[j], orbitals[k]);
					H.put(j, k, Huk);
					H.put(k, j, Huk);
				}
			}
		}
		
		//for (int j = 0; j < orbitals.length; j++) {
			//for (int k = 0; k < j + 1; k++) {
				//System.out.printf ("%-15f", H.get(j, k));
			//}
			//System.out.println ("");
		//}
		
		System.out.println ("1-electron matrix elements evaluated - moving on to two-electron matrix");

		//The idea of the integralarray is to simply store all the integrals in order they are called. It's basically my way of avoiding having to perform a Yoshemine sort.
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					
					for (int l: index[atomnumber[j]]) {
						if (l > -1) {
							MoleculeRun.
							integralarray.add(MNDO6G.getG(orbitals[j], orbitals[j], orbitals[l], orbitals[l]) - 0.5 * MNDO6G.getG(orbitals[j], orbitals[l], orbitals[j], orbitals[l]));
						}
					}
					
					for (int l: missingindex[atomnumber[j]]) {
						if (l > -1) {
							for (int m: missingindex[atomnumber[j]]) {
								if (m > -1) {
									if (atomnumber[l] == atomnumber[m]) {
										integralarray.add(MNDO6G.getG(orbitals[j], orbitals[j], orbitals[l], orbitals[m]));
									}
								}
								
							}
						}
					}
				}
				
				else if (atomnumber[j] == atomnumber[k]) {
					integralarray.add(1.5 * MNDO6G.getG(orbitals[j], orbitals[k], orbitals[j], orbitals[k]) - 0.5 * MNDO6G.getG(orbitals[j], orbitals[j], orbitals[k], orbitals[k]));
					
					for (int l: missingindex[atomnumber[j]]) {
						if (l > -1) {
							for (int m: missingindex[atomnumber[j]]) {
								if (m > -1) {
									if (atomnumber[l] == atomnumber[m]) {
										integralarray.add(MNDO6G.getG(orbitals[j], orbitals[k], orbitals[l], orbitals[m]));
									}
								}
								
							}
						}
					}
				}
				
				else {
					for (int l: index[atomnumber[j]]) {
						if (l > -1) {
							for (int m: index[atomnumber[k]]) {
								if (m > -1) {
									integralarray.add(-0.5 * MNDO6G.getG(orbitals[j], orbitals[l], orbitals[k], orbitals[m]));
								}
							}
						}
					}
				}
			}
		}
		System.out.println ("2-electron integrals evaluated");
		
		DoubleMatrix[] matrices = Eigen.symmetricEigenvectors(H);
		
		System.out.println ("initial diagonalization completed, beginning SCF iterations...");
		
		E = matrices[1].diag();
		
		C = matrices[0].transpose();
		
		G = DoubleMatrix.zeros (C.rows, C.columns);
		
		densitymatrix = DensityMatrix (C);
		
		DoubleMatrix olddensity = DoubleMatrix.zeros (C.rows, C.columns);
		
		F = H.dup();
		
		int integralcount;
		
		int numIt = 0;
		
		while (!isSimilar (densitymatrix, olddensity, 1E-10)) {//density matrix convergence criteria; since each iteration takes place within a fraction of a second I figured why not
			
			numIt++;
			olddensity = densitymatrix.dup();
			
			integralcount = 0;
			
			//this entire block of code fills up the G matrix, and it calls the integralarray to save time.
			
			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {
						
						for (int l: index[atomnumber[j]]) {
							if (l > -1) {
								val += densitymatrix.get(l,  l) * integralarray.get(integralcount);
								integralcount++;
							}
						}
						
						for (int l: missingindex[atomnumber[j]]) {
							if (l > -1) {
								for (int m: missingindex[atomnumber[j]]) {
									if (m > -1) {
										if (atomnumber[l] == atomnumber[m]) {
											val += densitymatrix.get(l,  m) * integralarray.get(integralcount);
											integralcount++;
										}
									}
									
								}
							}
						}
					}
					
					else if (atomnumber[j] == atomnumber[k]) {
						val += densitymatrix.get(j,  k) * integralarray.get(integralcount);
						integralcount++;
						
						for (int l: missingindex[atomnumber[j]]) {
							if (l > -1) {
								for (int m: missingindex[atomnumber[j]]) {
									if (m > -1) {
										if (atomnumber[l] == atomnumber[m]) {
											val += densitymatrix.get(l,  m) * integralarray.get(integralcount);
											integralcount++;
										}
									}
									
								}
							}
						}
					}
					
					else {
						for (int l: index[atomnumber[j]]) {
							if (l > -1) {
								for (int m: index[atomnumber[k]]) {
									if (m > -1) {
										val += densitymatrix.get(l,  m) * integralarray.get(integralcount);  
										integralcount++;
									}
								}
							}
						}
					}
					
					G.put(j,  k, val);
					G.put(k,  j, val);
				}
			}
			
			F = H.dup().add(G);
			
			
			matrices = Eigen.symmetricEigenvectors(F);
			
			E = matrices[1].diag();
			
			//System.err.println (E);
			
			C = matrices[0].transpose();
			
			densitymatrix = DensityMatrix (C).mmul(1-damp).add(olddensity.mmul(damp));
			
			//System.out.println (densitymatrix);
			
			if (numIt >= 10000) {
				System.err.println ("SCF Has Not Converged");
				
				System.err.println ("Damping Coefficient will be Increased, and the run restarted...");
				
				damp += 0.02;
				
				matrices = Eigen.symmetricEigenvectors(H);
				
				
				E = matrices[1].diag();
				
				C = matrices[0].transpose();
				
				G = DoubleMatrix.zeros (C.rows, C.columns);
				
				densitymatrix = DensityMatrix (C);
				
				numIt = 0;
				
				if (damp >= 1) {
					System.err.println ("Damping Coefficient Cannot Be Increased Further. Exiting program...");

					for (MNDOAtom a: atoms) {
						System.out.println (a.getZ() + "; " + Arrays.toString(a.getCoordinates()));
					}
					System.exit(0);
					
				}
			}
			
		}
		
		System.out.println ("SCF completed");
		
		//System.out.println ("Eigenvalues: " + E + "\n\n");
		
		//System.out.println ("Ionization energy: " + 0.001 * Math.round(E.get(NElectrons / 2 - 1, 0) * 1000) + " eV");
		
		//System.out.println ("LUMO energy: " + 0.001 * Math.round(E.get(NElectrons / 2, 0) * 1000) + " eV");
		
		double e = 0;
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * densitymatrix.get(j, k) * (H.get(j, k) + F.get(j, k));
			}
		}
		//System.out.println ("Electronic energy: " + 0.01 * Math.round(e * 100) + " eV");
		
		double d = 0;
		
		double heat = 0;
		
		for (int j = 0; j < atoms.length; j++) {
			heat += atoms[j].heat() - atoms[j].eisol();
			for (int k = j + 1; k < atoms.length; k++) {
				d += MNDOAtom.crf(atoms[j], atoms[k]);
				e += MNDOAtom.crf(atoms[j], atoms[k]);
			}
		}
		
		//System.out.println ("Core repulsion energy: " + 0.01 * Math.round(d * 100) + " eV");
		
		//System.out.println ("Energy: " + 0.01 * Math.round(e * 100) + " eV");
		
		energy = e;
		
		heat += e;
		
		this.hf =  heat / 4.3363E-2;
		
		if (NElectrons > 0) {
			this.homo = E.get(NElectrons / 2 - 1, 0);
		}
		else {
			this.homo = 0;
		}
		
		this.lumo = E.get(NElectrons / 2, 0);
		
		//System.out.println ("Heat of Formation: " + 0.01 * Math.round(heat * 100) + " eV = " + 0.01 * Math.round(heat * 100 / 4.3363E-2) + "kcal/mol");
		
		double[] populations  = new double[atoms.length];
		
		for (int j = 0; j < atoms.length; j++) {
			double sum = 0;
			for (int k: index[j]) {
				if (k > -1) {
					sum += densitymatrix.get(k, k);
				}
			}
			
			populations[j] = atoms[j].getQ() - sum;
		}
		
		
		double[] com = new double[] {0, 0, 0};
		
		double mass = 0;
		
		for (int j = 0; j < atoms.length; j++) {
			com[0] = com[0] + atoms[j].mass() * atoms[j].getCoordinates()[0];
			com[1] = com[1] + atoms[j].mass() * atoms[j].getCoordinates()[1];
			com[2] = com[2] + atoms[j].mass() * atoms[j].getCoordinates()[2];
			mass += atoms[j].mass();
		}
		
		com[0] = com[0] / mass;
		com[1] = com[1] / mass;
		com[2] = com[2] / mass;
		
		
		chargedip = new double[] {0, 0, 0};
		
		for (int j = 0; j < atoms.length; j++) {
			chargedip[0] += 2.5416 * populations[j] * (atoms[j].getCoordinates()[0] - com[0]); 
			chargedip[1] += 2.5416 * populations[j] * (atoms[j].getCoordinates()[1] - com[1]); 
			chargedip[2] += 2.5416 * populations[j] * (atoms[j].getCoordinates()[2] - com[2]); 
		}
		
		//chargedip[0] = 0.01 * Math.round(chargedip[0] * 100);
		//chargedip[1] = 0.01 * Math.round(chargedip[1] * 100);
		//chargedip[2] = 0.01 * Math.round(chargedip[2] * 100);
		
		//System.out.println ("point charge dipole contribution: " + Arrays.toString(chargedip));
		
		hybridip = new double[] {0, 0, 0};
		
		for (int j = 0; j < atoms.length; j++) {
			
			if (index[j][1] != -1) {//exclude hydrogen
				hybridip[0] = hybridip[0] -2.5416 * 2 * atoms[j].D1 * densitymatrix.get(index[j][0], index[j][1]);
				hybridip[1] = hybridip[1] -2.5416 * 2 * atoms[j].D1 * densitymatrix.get(index[j][0], index[j][2]);
				hybridip[2] = hybridip[2] -2.5416 * 2 * atoms[j].D1 * densitymatrix.get(index[j][0], index[j][3]);
			}
		}
		
		//hybridip[0] = 0.01 * Math.round(hybridip[0] * 100);
		//hybridip[1] = 0.01 * Math.round(hybridip[1] * 100);
		//hybridip[2] = 0.01 * Math.round(hybridip[2] * 100);
		
		//System.out.println ("hybrid dipole contribution: " + Arrays.toString(hybridip));
		
		dipoletot = new double[] {chargedip[0] + hybridip[0], chargedip[1] + hybridip[1], chargedip[2] + hybridip[2]};
		
		//System.out.println ("sum: " + Arrays.toString(dipoletot));
		
		
		dipole = Math.sqrt(dipoletot[0] * dipoletot[0] + dipoletot[1] * dipoletot[1] + dipoletot[2] * dipoletot[2]);
		
		//System.out.println ("dipole moment: " + dipole + " debyes");
		
	}
	
	public DoubleMatrix getE() {
		return E;
	}
	
	
	
	private DoubleMatrix DensityMatrix(DoubleMatrix c) {//density matrix construction by definition.
		
		
		DoubleMatrix densitymatrix = DoubleMatrix.zeros (orbitals.length, orbitals.length);
		
		for (int i = 0; i < orbitals.length; i++) {
			for (int j = 0; j < orbitals.length; j++) {
				double sum = 0;
				int count = NElectrons;
				int counter = -1;
				
				
				
				while (count > 0) {
					
					counter++;
					
					sum += 2 * c.get(counter, i) * c.get(counter, j);
					count -= 2;
					
				}
				
				
				densitymatrix.put(i,  j, sum);
			}
		}
		
		
		return densitymatrix;
	}
	
	private static boolean isSimilar (DoubleMatrix x, DoubleMatrix y, double limit) {
		
		for (int i = 0; i < y.rows; i++) {
			
			for (int j = 0; j < y.columns; j++) {
				
				if (Math.abs(x.get(i, j) - y.get(i, j)) > limit) {
					
					return false;
				}
			}
		}
		
		return true;
	}
	
	public DoubleMatrix densitymatrix() {
		return this.densitymatrix;
	}
}
