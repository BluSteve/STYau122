package nddoparam;

import org.ejml.simple.SimpleMatrix;
import org.jblas.DoubleMatrix;
import runcycle.input.RawMolecule;
import scf.Utils;

import java.util.Arrays;


public class SolutionU extends Solution {
	private SimpleMatrix Fa, Fb;
	private SimpleMatrix alphaDensity, betaDensity;

	protected SolutionU(NDDOAtom[] atoms, RawMolecule rm) {
		super(atoms, rm);
		this.mult = rm.mult;
		if (nElectrons % 2 == mult % 2 || mult < 1) {
			System.err.println(
					rm.debugName() + " Please check multiplicity and charge: " + nElectrons +
							", " + mult);
		}
		nElectrons -= (mult - 1);
	}

	@Override
	public SolutionU compute() {
		int nalpha = nElectrons / 2 + (mult - 1);
		int nbeta = nElectrons / 2;

		double[] integralArrayCoulomb = new double[rm.nCoulombInts];

		int integralCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							integralArrayCoulomb[integralCount] =
									NDDO6G.OneCenterERI(orbitals[j],
											orbitals[j],
											orbitals[l], orbitals[l]);
							integralCount++;
						}
					}
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (m > -1) {
									if (atomOfOrb[l] == atomOfOrb[m]) {
										integralArrayCoulomb[integralCount] =
												NDDO6G.getG(orbitals[j],
														orbitals[j],
														orbitals[l],
														orbitals[m]);
										integralCount++;
									}
								}
							}
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					integralArrayCoulomb[integralCount] = 2 *
							NDDO6G.OneCenterERI(orbitals[j], orbitals[k],
									orbitals[j],
									orbitals[k]);
					integralCount++;
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (m > -1) {
									if (atomOfOrb[l] == atomOfOrb[m]) {
										integralArrayCoulomb[integralCount] =
												NDDO6G.getG(orbitals[j],
														orbitals[k],
														orbitals[l],
														orbitals[m]);
										integralCount++;
									}
								}
							}
						}
					}
				}
			}
		}

		double[] integralArrayExchange = new double[rm.nExchangeInts];
		integralCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {

					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							integralArrayExchange[integralCount] = -1 *
									NDDO6G.OneCenterERI(orbitals[j],
											orbitals[l],
											orbitals[j], orbitals[l]);
							integralCount++;
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					integralArrayExchange[integralCount] = -1 *
							NDDO6G.OneCenterERI(orbitals[j], orbitals[k],
									orbitals[j],
									orbitals[k]) - 1 *
							NDDO6G.OneCenterERI(orbitals[j], orbitals[j],
									orbitals[k],
									orbitals[k]);
					integralCount++;
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
								if (m > -1) {
									integralArrayExchange[integralCount] = -1 *
											NDDO6G.getG(orbitals[j],
													orbitals[l],
													orbitals[k], orbitals[m]);
									integralCount++;
								}
							}
						}
					}
				}
			}
		}

		SimpleMatrix[] matrices = Utils.symEigen(H);

		SimpleMatrix ea = matrices[1].diag();

		SimpleMatrix eb = matrices[1].diag();

		SimpleMatrix ca = matrices[0].transpose();

		SimpleMatrix cb = matrices[0].transpose();

		SimpleMatrix j1 = new SimpleMatrix(ca.numRows(), ca.numCols());

		SimpleMatrix ka = new SimpleMatrix(ca.numRows(), ca.numCols());

		SimpleMatrix kb = new SimpleMatrix(ca.numRows(), ca.numCols());

		alphaDensity = calculateDensityMatrix(ca, nalpha);

		betaDensity = calculateDensityMatrix(cb, nbeta);

		SimpleMatrix oldalphadensity = Utils.filled(ca.numRows(), ca.numCols(),0);

		SimpleMatrix oldbetadensity = Utils.filled(ca.numRows(), ca.numCols(), 0);

		int Jcount, Kcount;

		int numIt = 0;
		boolean unstable = false;

		while (!(isSimilar(alphaDensity, oldalphadensity, 1E-10) &&
				isSimilar(betaDensity, oldbetadensity, 1E-10))) {

			numIt++;
			oldalphadensity = alphaDensity.copy();
			oldbetadensity = betaDensity.copy();

			Jcount = 0;
			Kcount = 0;

			//construct J matrix

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {

						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								val += (alphaDensity.get(l, l) +
										betaDensity.get(l, l)) *
										integralArrayCoulomb[Jcount];
								Jcount++;
							}
						}

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : missingOfAtom[atomOfOrb[j]]) {
									if (m > -1) {
										if (atomOfOrb[l] == atomOfOrb[m]) {
											val += (alphaDensity.get(l, m) +
													betaDensity.get(l, m)) *
													integralArrayCoulomb[Jcount];
											Jcount++;
										}
									}

								}
							}
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						val += (alphaDensity.get(j, k) +
								betaDensity.get(j, k)) *
								integralArrayCoulomb[Jcount];
						Jcount++;

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : missingOfAtom[atomOfOrb[j]]) {
									if (m > -1) {
										if (atomOfOrb[l] == atomOfOrb[m]) {
											val += (alphaDensity.get(l, m) +
													betaDensity.get(l, m)) *
													integralArrayCoulomb[Jcount];
											Jcount++;
										}
									}

								}
							}
						}
					}


					j1.set(j, k, val);
					j1.set(k, j, val);
				}
			}

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double vala = 0;
					double valb = 0;
					if (j == k) {

						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								vala += alphaDensity.get(l, l) *
										integralArrayExchange[Kcount];
								valb += betaDensity.get(l, l) *
										integralArrayExchange[Kcount];
								Kcount++;
							}
						}

					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						vala += alphaDensity.get(j, k) *
								integralArrayExchange[Kcount];
						valb += betaDensity.get(j, k) *
								integralArrayExchange[Kcount];
						Kcount++;

					}
					else {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : orbsOfAtom[atomOfOrb[k]]) {
									if (m > -1) {
										vala += alphaDensity.get(l, m) *
												integralArrayExchange[Kcount];
										valb += betaDensity.get(l, m) *
												integralArrayExchange[Kcount];
										Kcount++;
									}
								}
							}
						}
					}

					ka.set(j, k, vala);
					ka.set(k, j, vala);
					kb.set(j, k, valb);
					kb.set(k, j, valb);
				}
			}

			Fa = H.plus(j1).plus(ka);
			Fb = H.plus(j1).plus(kb);

			SimpleMatrix[] matrices1 = Utils.symEigen(Fa);

			SimpleMatrix[] matrices2 = Utils.symEigen(Fb);

			ea = matrices1[1].diag();

			eb = matrices2[1].diag();

			ca = matrices1[0].transpose();

			cb = matrices2[0].transpose();
			if (unstable) System.err.println(ea);

			double damp = 0.8;
			alphaDensity = calculateDensityMatrix(ca, nalpha).scale(1 - damp)
					.plus(oldalphadensity.scale(damp));

			betaDensity = calculateDensityMatrix(cb, nbeta).scale(1 - damp)
					.plus(oldbetadensity.scale(damp));

			if (numIt >= 1000000) {
				unstable = true;
				System.err.println("SCF Has Not Converged");

				System.err.println(
						"Damping Coefficient will be Increased, and the run " +
								"restarted." +
								"." +
								".");

				damp += 0.02;

				matrices = Utils.symEigen(H);

				System.out.println(
						"Exchange (K) matrix ERIs evaluated, beginning SCF " +
								"iterations." +
								"." +
								".");

				ea = matrices[1].diag();

				eb = matrices[1].diag();

				ca = matrices[0].transpose();

				cb = matrices[0].transpose();

				j1 = new SimpleMatrix(ca.numRows(), ca.numCols());

				ka = new SimpleMatrix(ca.numRows(), ca.numCols());

				kb = new SimpleMatrix(ca.numRows(), ca.numCols());

				alphaDensity = calculateDensityMatrix(ca, nalpha);

				betaDensity = calculateDensityMatrix(cb, nbeta);

				numIt = 0;

				if (damp >= 1) {
					System.err.println(
							"Damping Coefficient Cannot Be Increased Further" +
									"." +
									" " +
									"Exiting " +
									"program...");

					for (NDDOAtom a : atoms) {
						System.out.println(a.getAtomProperties().getZ() + ";" +
								" " +
								Arrays.toString(a.getCoordinates()));
					}

				}
			}

		}


		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * alphaDensity.get(j, k) *
						(H.get(j, k) + Fa.get(j, k));
				e += 0.5 * betaDensity.get(j, k) * (H.get(j, k) + Fb.get(j,
						k));
			}
		}
		double checksum = 0;

		for (int a = 0; a < atoms.length; a++) {
			checksum += E(a, orbsOfAtom);
		}

		for (int a = 0; a < atoms.length; a++) {
			for (int b = a + 1; b < atoms.length; b++) {
				checksum += E(a, b, orbsOfAtom);
			}
		}

		if (Math.abs(checksum - e) > 1E-5 || checksum != checksum) {
			System.err.println("I knew it!");
			System.err.println(checksum);
			System.err.println(e);
//			System.exit(0);
		}


		double heat = 0;
		for (int j = 0; j < atoms.length; j++) {
			heat += atoms[j].getHeat() - atoms[j].getParams().getEisol();
			for (int k = j + 1; k < atoms.length; k++) {
				e += atoms[j].crf(atoms[k]);
			}
		}
		energy = e;

		heat += e;

		this.hf = heat / 4.3363E-2;

		this.homo = ea.get(nalpha - 1, 0);
		this.lumo = 0.001 * Math.round(eb.get(nbeta, 0) * 1000);
//		System.out.println(moleculeName + " SCF completed");

		double[] populations = new double[atoms.length];

		for (int j = 0; j < atoms.length; j++) {
			double sum = 0;
			for (int k : orbsOfAtom[j]) {
				if (k > -1) {
					sum += alphaDensity.get(k, k) + betaDensity.get(k, k);
				}
			}

			populations[j] = atoms[j].getAtomProperties().getQ() - sum;
		}


		double[] com = new double[]{0, 0, 0};

		double mass = 0;

		for (NDDOAtom atom : atoms) {
			com[0] = com[0] + atom.getMass() * atom.getCoordinates()[0];
			com[1] = com[1] + atom.getMass() * atom.getCoordinates()[1];
			com[2] = com[2] + atom.getMass() * atom.getCoordinates()[2];
			mass += atom.getMass();
		}

		com[0] = com[0] / mass;
		com[1] = com[1] / mass;
		com[2] = com[2] / mass;


		chargedip = new double[]{0, 0, 0};

		for (int j = 0; j < atoms.length; j++) {
			chargedip[0] +=
					2.5416 * populations[j] *
							(atoms[j].getCoordinates()[0] - com[0]);
			chargedip[1] +=
					2.5416 * populations[j] *
							(atoms[j].getCoordinates()[1] - com[1]);
			chargedip[2] +=
					2.5416 * populations[j] *
							(atoms[j].getCoordinates()[2] - com[2]);
		}


		hybridip = new double[]{0, 0, 0};

		for (int j = 0; j < atoms.length; j++) {

			if (orbsOfAtom[j][1] != -1) {//exclude hydrogen
				hybridip[0] = hybridip[0] - 2.5416 * 2 * atoms[j].D1 *
						(alphaDensity.get(orbsOfAtom[j][0], orbsOfAtom[j][1]) +
								betaDensity.get(orbsOfAtom[j][0],
										orbsOfAtom[j][1]));
				hybridip[1] = hybridip[1] - 2.5416 * 2 * atoms[j].D1 *
						(alphaDensity.get(orbsOfAtom[j][0], orbsOfAtom[j][2]) +
								betaDensity.get(orbsOfAtom[j][0],
										orbsOfAtom[j][2]));
				hybridip[2] = hybridip[2] - 2.5416 * 2 * atoms[j].D1 *
						(alphaDensity.get(orbsOfAtom[j][0], orbsOfAtom[j][3]) +
								betaDensity.get(orbsOfAtom[j][0],
										orbsOfAtom[j][3]));
			}
		}


		dipoletot = new double[]{chargedip[0] + hybridip[0],
				chargedip[1] + hybridip[1],
				chargedip[2] + hybridip[2]};


		dipole = Math.sqrt(
				dipoletot[0] * dipoletot[0] + dipoletot[1] * dipoletot[1] +
						dipoletot[2] * dipoletot[2]);
		return this;
	}

	@SuppressWarnings("DuplicatedCode")
	protected int findNCoulombInts() {
		int size = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							size++;
						}
					}
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (m > -1) {
									if (atomOfOrb[l] == atomOfOrb[m]) {
										size++;
									}
								}
							}
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					size++;
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (m > -1) {
									if (atomOfOrb[l] == atomOfOrb[m]) {
										size++;
									}
								}

							}
						}
					}
				}
			}
		}
		return size;
	}

	protected int findNExchangeInts() {
		int size = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							size++;
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					size++;
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
								if (m > -1) {
									size++;
								}
							}
						}
					}
				}
			}
		}
		return size;
	}

	private double E(int atomnum, int[][] index) {

		SimpleMatrix densitymatrix = alphaDensity.plus(betaDensity);

		double e = 0;

		for (int i : index[atomnum]) {
			if (i > -1) {
				e += densitymatrix.get(i, i) * orbitals[i].U();
			}
		}

		for (int i : index[atomnum]) {
			for (int j : index[atomnum]) {
				for (int k : index[atomnum]) {
					for (int l : index[atomnum]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {

							e += 0.5 * densitymatrix.get(i, j) *
									densitymatrix.get(k,
											l) *
									NDDO6G.OneCenterERI(orbitals[i],
											orbitals[j],
											orbitals[k], orbitals[l]);
							e -= 0.5 * (alphaDensity.get(i, j) *
									alphaDensity.get(k, l) +
									betaDensity.get(i, j) *
											betaDensity.get(k, l)) *
									NDDO6G.OneCenterERI(orbitals[i],
											orbitals[k],
											orbitals[j], orbitals[l]);
						}
					}
				}
			}
		}

		return e;
	}


	private double E(int atomnum1, int atomnum2, int[][] index) {

		SimpleMatrix densitymatrix = alphaDensity.plus(betaDensity);

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += densitymatrix.get(i, j) *
							atoms[atomnum2].V(orbitals[i], orbitals[j]);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e += densitymatrix.get(k, l) *
							atoms[atomnum1].V(orbitals[k], orbitals[l]);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * densitymatrix.get(i, k) *
							NDDO6G.beta(orbitals[i], orbitals[k]);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {

							e += (densitymatrix.get(i, j) *
									densitymatrix.get(k, l) -
									alphaDensity.get(i, k) *
											alphaDensity.get(j, l) -
									betaDensity.get(i, k) *
											betaDensity.get(j, l))
									* NDDO6G.getG(orbitals[i], orbitals[j],
									orbitals[k],
									orbitals[l]);
						}
					}
				}
			}
		}

		return e;


	}

	private SimpleMatrix calculateDensityMatrix(SimpleMatrix c,
												int NElectrons) {
		SimpleMatrix densityMatrix =
				new SimpleMatrix(orbitals.length, orbitals.length);

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
				densityMatrix.set(i, j, sum);
			}
		}
		return densityMatrix;
	}

	@Override
	public DoubleMatrix alphaDensity() {
		return Utils.toDoubleMatrix(alphaDensity);
	}

	@Override
	public DoubleMatrix betaDensity() {
		return Utils.toDoubleMatrix(betaDensity);
	}

	@Override
	public DoubleMatrix densityMatrix() {
		return Utils.toDoubleMatrix(alphaDensity.plus(betaDensity));
	}
}
