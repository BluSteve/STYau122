package nddo.solution;

import nddo.NDDO6G;
import nddo.NDDOAtom;
import org.ejml.simple.SimpleMatrix;
import runcycle.input.RawMolecule;
import tools.Utils;


public class SolutionU extends Solution {
	private SimpleMatrix Fa, Fb, alphaDensity, betaDensity;

	protected SolutionU(NDDOAtom[] atoms, RawMolecule rm) {
		super(atoms, rm);
		if (nElectrons % 2 == mult % 2 || mult < 1) {
			rm.getLogger().error("Please check mult and charge: " +
					"nElectrons: {}, mult: {}", nElectrons, mult);
		}
		nElectrons -= mult - 1;
	}

	@Override
	public SolutionU compute() {
		double damp = 0.8;
		int nalpha = nElectrons / 2 + (mult - 1);
		int nbeta = nElectrons / 2;

		double[] integralArrayCoulomb = new double[getRm().nCoulombInts];
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

		double[] integralArrayExchange = new double[getRm().nExchangeInts];
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

		SimpleMatrix oldalphadensity =
				new SimpleMatrix(ca.numRows(), ca.numCols());
		SimpleMatrix oldbetadensity =
				new SimpleMatrix(ca.numRows(), ca.numCols());

		int Jcount, Kcount;
		int numIt = 0;

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

			SimpleMatrix[] matricesa = Utils.symEigen(Fa);
			SimpleMatrix[] matricesb = Utils.symEigen(Fb);

			ea = matricesa[1].diag();
			eb = matricesb[1].diag();
			ca = matricesa[0].transpose();
			cb = matricesb[0].transpose();

			alphaDensity = calculateDensityMatrix(ca, nalpha).scale(1 - damp)
					.plus(oldalphadensity.scale(damp));
			betaDensity = calculateDensityMatrix(cb, nbeta).scale(1 - damp)
					.plus(oldbetadensity.scale(damp));

			if (numIt >= 1000000) {
				getRm().getLogger().error("SCF has not converged. " +
						"Damping coefficient (currently {}) will be " +
						"increased, and the run restarted.", damp);

				damp += 0.02;

				matrices = Utils.symEigen(H);

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
					getRm().getLogger().error(
							"Damping coefficient cannot be increased" +
									" further.");
				}
			}
		}

		findEnergyAndHf();
		this.homo = ea.get(nalpha - 1, 0);
		this.lumo = 0.001 * Math.round(eb.get(nbeta, 0) * 1000);
		findDipole();

		return this;
	}

	private void findEnergyAndHf() {
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				energy += 0.5 * alphaDensity.get(j, k) *
						(H.get(j, k) + Fa.get(j, k));
				energy += 0.5 * betaDensity.get(j, k) *
						(H.get(j, k) + Fb.get(j, k));
			}
		}

		double heat = 0;
		for (int j = 0; j < atoms.length; j++) {
			heat += atoms[j].getHeat() - atoms[j].getParams().getEisol();
			for (int k = j + 1; k < atoms.length; k++) {
				energy += atoms[j].crf(atoms[k]);
			}
		}

		heat += energy;
		hf = heat / 4.3363E-2;
	}

	private void findDipole() {
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
			com[0] += atom.getMass() * atom.getCoordinates()[0];
			com[1] += atom.getMass() * atom.getCoordinates()[1];
			com[2] += atom.getMass() * atom.getCoordinates()[2];
			mass += atom.getMass();
		}

		com[0] /= mass;
		com[1] /= mass;
		com[2] /= mass;

		chargedip = new double[]{0, 0, 0};

		double v = 2.5416;
		for (int j = 0; j < atoms.length; j++) {
			double v1 = v * populations[j];
			chargedip[0] += v1 *
					(atoms[j].getCoordinates()[0] - com[0]);
			chargedip[1] += v1 *
					(atoms[j].getCoordinates()[1] - com[1]);
			chargedip[2] += v1 *
					(atoms[j].getCoordinates()[2] - com[2]);
		}

		hybridip = new double[]{0, 0, 0};

		for (int j = 0; j < atoms.length; j++) {
			if (orbsOfAtom[j][1] != -1) { // exclude hydrogen
				double v1 = v * 2 * atoms[j].D1;
				hybridip[0] -= v1 *
						(alphaDensity.get(orbsOfAtom[j][0], orbsOfAtom[j][1]) +
								betaDensity.get(orbsOfAtom[j][0],
										orbsOfAtom[j][1]));
				hybridip[1] -= v1 *
						(alphaDensity.get(orbsOfAtom[j][0], orbsOfAtom[j][2]) +
								betaDensity.get(orbsOfAtom[j][0],
										orbsOfAtom[j][2]));
				hybridip[2] -= v1 *
						(alphaDensity.get(orbsOfAtom[j][0], orbsOfAtom[j][3]) +
								betaDensity.get(orbsOfAtom[j][0],
										orbsOfAtom[j][3]));
			}
		}

		dipoletot = new double[]{chargedip[0] + hybridip[0],
				chargedip[1] + hybridip[1],
				chargedip[2] + hybridip[2]};

		dipole = Math.sqrt(dipoletot[0] * dipoletot[0] +
				dipoletot[1] * dipoletot[1] +
				dipoletot[2] * dipoletot[2]);
	}

	private SimpleMatrix calculateDensityMatrix(SimpleMatrix c,
												int nElectrons) {
		SimpleMatrix densityMatrix = new SimpleMatrix(orbitals.length,
				orbitals.length);

		for (int i = 0; i < orbitals.length; i++) {
			for (int j = 0; j < orbitals.length; j++) {
				double sum = 0;
				int count = nElectrons;
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

	@Override
	public SimpleMatrix alphaDensity() {
		return alphaDensity;
	}

	@Override
	public SimpleMatrix betaDensity() {
		return betaDensity;
	}

	@Override
	public SimpleMatrix densityMatrix() {
		return alphaDensity.plus(betaDensity);
	}
}
