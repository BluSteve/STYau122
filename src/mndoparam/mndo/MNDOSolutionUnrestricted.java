package mndoparam.mndo;

import java.util.ArrayList;
import java.util.Arrays;

import org.jblas.DoubleMatrix;
import org.jblas.Eigen;


public class MNDOSolutionUnrestricted {
	private DoubleMatrix H, J, Ka, Kb, Fa, Fb, Ca, Cb;

	DoubleMatrix Ea;

	DoubleMatrix Eb;

	private DoubleMatrix alphadensity;

	private DoubleMatrix betadensity;
	
	private MNDOAtom[] atoms;
	
	private MNDO6G[] orbitals;
	
	private int Nalpha, Nbeta;
	
	private ArrayList<Double> integralarraycoulomb, integralarrayexchange;
	
	public double energy, homo, lumo, hf;
	
	public double[] chargedip, hybridip, dipoletot;
	
	public double dipole;
	
	public double electronicenergy;
	
	public int charge, multiplicity;
	
	public MNDOSolutionUnrestricted (MNDOAtom[] atoms, int charge, int multiplicity) {
		
		this.charge = charge;
		
		this.multiplicity = multiplicity;
		
		double damp = 0.8;
		
		int n = 0;
		
		for (MNDOAtom a: atoms) {
			n += a.getQ();
		}
		
		n -= charge;
		
		if (n % 2 == multiplicity % 2 || multiplicity < 1) {
			System.err.println ("You're high. (Please check multiplicity and charge)");
			System.exit(0);
		}
		
		n -= (multiplicity - 1);
		
		Nalpha = n / 2 + (multiplicity - 1);
		
		Nbeta = n / 2;
		
		
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
				//System.err.println ("(" + j + ", " + k + ")");
				if (k == j) {
					double Huu = orbitals[j].U();
					//System.err.println ("U: " + Huu);
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
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < j + 1; k++) {
				//System.out.printf ("%-15f", H.get(j, k));
			}
			//System.out.println ("");
		}
		
		System.out.println ("1-electron matrix elements evaluated - moving on to two-electron matrix");
		
		integralarraycoulomb = new ArrayList<Double>();
		
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				//System.err.println ("(" + j + ", " + k + ")");
				if (j == k) {
					
					for (int l: index[atomnumber[j]]) {
						if (l > -1) {
							integralarraycoulomb.add(MNDO6G.OneCenterERI(orbitals[j], orbitals[j], orbitals[l], orbitals[l]));
						}
					}
					
					for (int l: missingindex[atomnumber[j]]) {
						if (l > -1) {
							for (int m: missingindex[atomnumber[j]]) {
								if (m > -1) {
									if (atomnumber[l] == atomnumber[m]) {
										integralarraycoulomb.add(MNDO6G.getG(orbitals[j], orbitals[j], orbitals[l], orbitals[m]));
									}
								}
								
							}
						}
					}
				}
				
				else if (atomnumber[j] == atomnumber[k]) {
					integralarraycoulomb.add(2 * MNDO6G.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]));
					
					for (int l: missingindex[atomnumber[j]]) {
						if (l > -1) {
							for (int m: missingindex[atomnumber[j]]) {
								if (m > -1) {
									if (atomnumber[l] == atomnumber[m]) {
										integralarraycoulomb.add(MNDO6G.getG(orbitals[j], orbitals[k], orbitals[l], orbitals[m]));
									}
								}
								
							}
						}
					}
				}
			}
		}
		
		System.out.println ("Coulomb (J) matrix ERIs evaluated - moving on to Exchange (K) matrix ERIs...");
		
		integralarrayexchange = new ArrayList<Double>();
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				//System.err.println ("(" + j + ", " + k + ")");
				if (j == k) {
					
					for (int l: index[atomnumber[j]]) {
						if (l > -1) {
							integralarrayexchange.add(- 1 * MNDO6G.OneCenterERI(orbitals[j], orbitals[l], orbitals[j], orbitals[l]));
						}
					}
				}
				
				else if (atomnumber[j] == atomnumber[k]) {
					//System.err.println ("1.5[" + j + k + "|" + j + k + "] - 0.5[" + j + j + "|" + k + k + "]");
					integralarrayexchange.add(-1 * MNDO6G.OneCenterERI(orbitals[j], orbitals[k], orbitals[j], orbitals[k]) - 1 * MNDO6G.OneCenterERI(orbitals[j], orbitals[j], orbitals[k], orbitals[k]));
				}
				
				else {
					for (int l: index[atomnumber[j]]) {
						if (l > -1) {
							for (int m: index[atomnumber[k]]) {
								if (m > -1) {
									integralarrayexchange.add(-1 * MNDO6G.getG(orbitals[j], orbitals[l], orbitals[k], orbitals[m]));
								}
							}
						}
					}
				}
			}
		}
		
		DoubleMatrix[] matrices = Eigen.symmetricEigenvectors(H);
		
		System.out.println ("Exchange (K) matrix ERIs evaluated, beginning SCF iterations...");
		
		Ea = matrices[1].diag();
		
		Eb = matrices[1].diag();
		
		Ca = matrices[0].transpose();
		
		Cb = matrices[0].transpose();
		
		J = new DoubleMatrix (Ca.rows, Ca.columns);
		
		Ka = new DoubleMatrix (Ca.rows, Ca.columns);
		
		Kb = new DoubleMatrix (Ca.rows, Ca.columns);
		
		alphadensity = DensityMatrix (Ca, Nalpha);
		
		betadensity = DensityMatrix (Cb, Nbeta);
		
		DoubleMatrix oldalphadensity = DoubleMatrix.zeros (Ca.rows, Ca.columns);
		
		DoubleMatrix oldbetadensity = DoubleMatrix.zeros (Ca.rows, Ca.columns);
		
		int Jcount = 0;
		
		int Kcount = 0;
		
		int numIt = 0;
		while (!(isSimilar (alphadensity, oldalphadensity, 1E-10) && isSimilar (betadensity, oldbetadensity, 1E-10))) {
			
			numIt++;
			oldalphadensity = alphadensity.dup();
			oldbetadensity = betadensity.dup();
			
			Jcount = 0;
			Kcount = 0;
			
			//construct J matrix
			
			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {
						
						for (int l: index[atomnumber[j]]) {
							if (l > -1) {
								val += (alphadensity.get(l,  l) + betadensity.get(l, l)) * integralarraycoulomb.get(Jcount);
								Jcount++;
							}
						}
						
						for (int l: missingindex[atomnumber[j]]) {
							if (l > -1) {
								for (int m: missingindex[atomnumber[j]]) {
									if (m > -1) {
										if (atomnumber[l] == atomnumber[m]) {
											val += (alphadensity.get(l,  m) + betadensity.get(l, m)) * integralarraycoulomb.get(Jcount);
											Jcount++;
										}
									}
									
								}
							}
						}
					}
					
					else if (atomnumber[j] == atomnumber[k]) {
						val += (alphadensity.get(j,  k) + betadensity.get(j, k)) * integralarraycoulomb.get(Jcount);
						Jcount++;
						
						for (int l: missingindex[atomnumber[j]]) {
							if (l > -1) {
								for (int m: missingindex[atomnumber[j]]) {
									if (m > -1) {
										if (atomnumber[l] == atomnumber[m]) {
											val += (alphadensity.get(l,  m) + betadensity.get(l, m)) * integralarraycoulomb.get(Jcount);
											Jcount++;
										}
									}
									
								}
							}
						}
					}

					
					J.put(j,  k, val);
					J.put(k,  j, val);
				}
			}
			
			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double vala = 0;
					double valb = 0;
					if (j == k) {
						
						for (int l: index[atomnumber[j]]) {
							if (l > -1) {
								vala += alphadensity.get(l,  l) * integralarrayexchange.get(Kcount);
								valb += betadensity.get(l,  l) * integralarrayexchange.get(Kcount);
								Kcount++;
							}
						}
						
					}
					
					else if (atomnumber[j] == atomnumber[k]) {
						vala += alphadensity.get(j,  k) * integralarrayexchange.get(Kcount);
						valb += betadensity.get(j,  k) * integralarrayexchange.get(Kcount);
						Kcount++;

					}
					
					else {
						for (int l: index[atomnumber[j]]) {
							if (l > -1) {
								for (int m: index[atomnumber[k]]) {
									if (m > -1) {
										vala += alphadensity.get(l,  m) * integralarrayexchange.get(Kcount);  
										valb += betadensity.get(l,  m) * integralarrayexchange.get(Kcount);  
										Kcount++;
									}
								}
							}
						}
					}
					
					Ka.put(j,  k, vala);
					Ka.put(k,  j, vala);
					Kb.put(j,  k, valb);
					Kb.put(k,  j, valb);
				}
			}
			
			Fa = H.add(J).add(Ka);
			Fb = H.add(J).add(Kb);
			
			DoubleMatrix[] matrices1 = Eigen.symmetricEigenvectors(Fa);
			
			DoubleMatrix[] matrices2 = Eigen.symmetricEigenvectors(Fb);
			
			Ea = matrices1[1].diag();
			
			Eb = matrices2[1].diag();
			
			Ca = matrices1[0].transpose();
			
			Cb = matrices2[0].transpose();
			
			alphadensity = DensityMatrix (Ca, Nalpha).mmul(1-damp).add(oldalphadensity.mmul(damp));
			
			betadensity = DensityMatrix (Cb, Nbeta).mmul(1-damp).add(oldbetadensity.mmul(damp));
			
			if (numIt >= 100000) {
				System.err.println ("SCF Has Not Converged");
				
				System.err.println ("Damping Coefficient will be Increased, and the run restarted...");
				
				damp += 0.02;
				
				matrices = Eigen.symmetricEigenvectors(H);
				
				matrices = Eigen.symmetricEigenvectors(H);
				
				System.out.println ("Exchange (K) matrix ERIs evaluated, beginning SCF iterations...");
				
				Ea = matrices[1].diag();
				
				Eb = matrices[1].diag();
				
				Ca = matrices[0].transpose();
				
				Cb = matrices[0].transpose();
				
				J = new DoubleMatrix (Ca.rows, Ca.columns);
				
				Ka = new DoubleMatrix (Ca.rows, Ca.columns);
				
				Kb = new DoubleMatrix (Ca.rows, Ca.columns);
				
				alphadensity = DensityMatrix (Ca, Nalpha);
				
				betadensity = DensityMatrix (Cb, Nbeta);
				
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
		
		//System.out.println ("Alpha eigenvalues: " + Ea);
		//System.out.println ("Beta eigenvalues: " + Eb);
		
		double e = 0;
		
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * alphadensity.get(j, k) * (H.get(j, k) + Fa.get(j, k));
				e += 0.5 * betadensity.get(j, k) * (H.get(j, k) + Fb.get(j, k));
			}
		}
		
		electronicenergy = e;
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
		
		this.hf =  heat/ 4.3363E-2;
		
		this.homo = Ea.get(Nalpha - 1, 0);
		this.lumo = 0.001 * Math.round(Eb.get(Nbeta, 0) * 1000);
		
		//System.out.println ("HOMO energy: " + homo + " eV");
		
		//System.out.println ("LUMO energy: " + lumo + " eV");
		
		//System.out.println ("Heat of Formation: " + 0.01 * Math.round(heat * 100) + " eV = " + 0.01 * Math.round(heat * 100 / 4.3363E-2) + "kcal/mol");
		
		double[] populations  = new double[atoms.length];
		
		for (int j = 0; j < atoms.length; j++) {
			double sum = 0;
			for (int k: index[j]) {
				if (k > -1) {
					sum += alphadensity.get(k, k) + betadensity.get(k, k);
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
		
		//chargedip[0] = 0.001 * Math.round(chargedip[0] * 1000);
		//chargedip[1] = 0.001 * Math.round(chargedip[1] * 1000);
		//chargedip[2] = 0.001 * Math.round(chargedip[2] * 1000);
		
		//System.out.println ("point charge dipole contribution: " + Arrays.toString(chargedip));
		
		hybridip = new double[] {0, 0, 0};
		
		for (int j = 0; j < atoms.length; j++) {
			
			if (index[j][1] != -1) {//exclude hydrogen
				hybridip[0] = hybridip[0] -2.5416 * 2 * atoms[j].D1 * (alphadensity.get(index[j][0], index[j][1]) + betadensity.get(index[j][0], index[j][1]));
				hybridip[1] = hybridip[1] -2.5416 * 2 * atoms[j].D1 * (alphadensity.get(index[j][0], index[j][2]) + betadensity.get(index[j][0], index[j][2]));
				hybridip[2] = hybridip[2] -2.5416 * 2 * atoms[j].D1 * (alphadensity.get(index[j][0], index[j][3]) + betadensity.get(index[j][0], index[j][3]));
			}
		}
		
		//hybridip[0] = 0.001 * Math.round(hybridip[0] * 1000);
		//hybridip[1] = 0.001 * Math.round(hybridip[1] * 1000);
		//hybridip[2] = 0.001 * Math.round(hybridip[2] * 1000);
		
		//System.out.println ("hybrid dipole contribution: " + Arrays.toString(hybridip));
		
		dipoletot = new double[] {chargedip[0] + hybridip[0], chargedip[1] + hybridip[1], chargedip[2] + hybridip[2]};
		
		//System.out.println ("sum: " + Arrays.toString(dipoletot));
		
		
		dipole = Math.sqrt(dipoletot[0] * dipoletot[0] + dipoletot[1] * dipoletot[1] + dipoletot[2] * dipoletot[2]);
		
		//System.out.println ("dipole moment: " + dipole + " debyes");
	}
	
	private DoubleMatrix DensityMatrix(DoubleMatrix c, int NElectrons) {//density matrix construction by definition.
		
		
		DoubleMatrix densitymatrix = new DoubleMatrix (orbitals.length, orbitals.length);
		
		for (int i = 0; i < orbitals.length; i++) {
			for (int j = 0; j < orbitals.length; j++) {
				double sum = 0;
				int count = NElectrons;
				int counter = -1;
				
				
				
				while (count > 0) {
					
					counter++;
					
					sum += c.get(counter, i) * c.get(counter, j);
					count -= 1;
					
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
	
	public DoubleMatrix alphadensity() {
		return this.alphadensity;
	}
	
	public DoubleMatrix betadensity() {
		return this.betadensity;
	}
}
