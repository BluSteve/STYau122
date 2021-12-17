package nddo.solution;

import nddo.NDDOAtom;
import nddo.structs.MoleculeInfo;
import org.ejml.simple.SimpleMatrix;
import tools.Utils;

import java.util.Arrays;

import static nddo.State.nom;


public class SolutionU extends Solution {
	public SimpleMatrix Ea, Eb, ca, cb, Fa;
	public double[] integralArrayCoulomb, integralArrayExchange;
	private SimpleMatrix Fb;
	private SimpleMatrix alphaDensity, betaDensity;

	protected SolutionU(MoleculeInfo mi, NDDOAtom[] atoms) {
		super(mi, atoms);
		if (nElectrons % 2 == mult % 2 || mult < 1) {
			mi.getLogger().error("Please check mult and charge: nElectrons: {}, mult: {}", nElectrons, mult);
		}
	}

	private static SimpleMatrix commutator(SimpleMatrix F, SimpleMatrix D) {
		return F.mult(D).minus(D.mult(F));
	}

	@Override
	public SolutionU compute() {
		int nalpha = (nElectrons - mult + 1) / 2 + mult - 1;
		int nbeta = nElectrons - nalpha;

		integralArrayCoulomb = new double[rm.nCoulombInts];

		int integralCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							integralArrayCoulomb[integralCount] =
									nom.OneCenterERI(orbitals[j],
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
												nom.getG(orbitals[j],
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
							nom.OneCenterERI(orbitals[j], orbitals[k],
									orbitals[j],
									orbitals[k]);
					integralCount++;
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (m > -1) {
									if (atomOfOrb[l] == atomOfOrb[m]) {
										integralArrayCoulomb[integralCount] =
												nom.getG(orbitals[j],
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

		integralArrayExchange = new double[rm.nExchangeInts];
		integralCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {

					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							integralArrayExchange[integralCount] = -1 *
									nom.OneCenterERI(orbitals[j],
											orbitals[l],
											orbitals[j], orbitals[l]);
							integralCount++;
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					integralArrayExchange[integralCount] = -1 *
							nom.OneCenterERI(orbitals[j], orbitals[k],
									orbitals[j],
									orbitals[k]) - 1 *
							nom.OneCenterERI(orbitals[j], orbitals[j],
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
											nom.getG(orbitals[j],
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

		Ea = matrices[1].diag();

		Eb = matrices[1].diag();

		ca = matrices[0].transpose();

		cb = matrices[0].transpose();

		SimpleMatrix j1 = new SimpleMatrix(ca.numRows(), ca.numCols());

		SimpleMatrix ka = new SimpleMatrix(ca.numRows(), ca.numCols());

		SimpleMatrix kb = new SimpleMatrix(ca.numRows(), ca.numCols());

		alphaDensity = calculateDensityMatrix(ca, nalpha);

		betaDensity = calculateDensityMatrix(cb, nbeta);

		SimpleMatrix oldalphadensity = new SimpleMatrix(ca.numRows(), ca.numCols());

		SimpleMatrix oldbetadensity = new SimpleMatrix(ca.numRows(), ca.numCols());

		int Jcount, Kcount;

		int numIt = 0;

		double damp = 0.55;

		SimpleMatrix[] Farrayalpha = new SimpleMatrix[8];
		SimpleMatrix[] Farraybeta = new SimpleMatrix[8];

		SimpleMatrix[] Darrayalpha = new SimpleMatrix[8];
		SimpleMatrix[] Darraybeta = new SimpleMatrix[8];


		SimpleMatrix B = new SimpleMatrix(8, 8);

		SimpleMatrix[] commutatorarrayalpha = new SimpleMatrix[8];
		SimpleMatrix[] commutatorarraybeta = new SimpleMatrix[8];

		double DIISError = 10;

		while (DIISError > 2 * 1E-13) {

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

			if (numIt < Farrayalpha.length) {

				Farrayalpha[numIt] = Fa.copy();
				Farraybeta[numIt] = Fb.copy();


				Darrayalpha[numIt] = alphaDensity.copy();
				Darraybeta[numIt] = betaDensity.copy();

				commutatorarrayalpha[numIt] = commutator(Fa.copy(), alphaDensity.copy());
				commutatorarraybeta[numIt] = commutator(Fb.copy(), betaDensity.copy());

				DIISError = commutatorarrayalpha[numIt].normF() + commutatorarraybeta[numIt].normF();

				for (int i = 0; i <= numIt; i++) {

					double product = commutatorarrayalpha[numIt].mult(commutatorarrayalpha[i].transpose()).diag()
									.elementSum() +
									commutatorarraybeta[numIt].mult(commutatorarraybeta[i].transpose()).diag()
											.elementSum();
					B.set(i, numIt, product);
					B.set(numIt, i, product);

				}
			}

			else {

				for (int i = 0; i < Farrayalpha.length - 1; i++) {

					Farrayalpha[i] = Farrayalpha[i + 1].copy();
					Farraybeta[i] = Farraybeta[i + 1].copy();

					Darrayalpha[i] = Darrayalpha[i + 1].copy();
					Darraybeta[i] = Darraybeta[i + 1].copy();

					commutatorarrayalpha[i] = commutatorarrayalpha[i + 1].copy();
					commutatorarraybeta[i] = commutatorarraybeta[i + 1].copy();
				}

				Farrayalpha[Farrayalpha.length - 1] = Fa.copy();
				Farraybeta[Farraybeta.length - 1] = Fb.copy();

				Darrayalpha[Darrayalpha.length - 1] = alphaDensity.copy();
				Darraybeta[Darraybeta.length - 1] = betaDensity.copy();

				commutatorarrayalpha[Darrayalpha.length - 1] = commutator(Fa.copy(), alphaDensity.copy());
				commutatorarraybeta[Darraybeta.length - 1] = commutator(Fb.copy(), betaDensity.copy());

				DIISError = commutatorarrayalpha[Darrayalpha.length - 1].normF() +
						commutatorarraybeta[Darraybeta.length - 1].normF();

				// B is dy/dx sort of, make dy/dx 0
				SimpleMatrix newB = new SimpleMatrix(8, 8);

				for (int i = 0; i < Farrayalpha.length - 1; i++) {
					for (int j = i; j < Farrayalpha.length - 1; j++) {
						newB.set(i, j, B.get(i + 1, j + 1));
						newB.set(j, i, B.get(i + 1, j + 1));
					}
				}

				for (int i = 0; i < Farrayalpha.length; i++) {

					double product =
							commutatorarrayalpha[Farrayalpha.length - 1].transpose().mult(commutatorarrayalpha[i])
									.diag().elementSum() +
									commutatorarraybeta[Farraybeta.length - 1].transpose().mult(commutatorarraybeta[i])
											.diag().elementSum();

					newB.set(i, Farrayalpha.length - 1, product);
					newB.set(Farrayalpha.length - 1, i, product);
				}

				B = newB.copy();
			}

			SimpleMatrix mat = new SimpleMatrix(Math.min(Farrayalpha.length + 1, numIt + 2),
					Math.min(Farrayalpha.length + 1, numIt + 2));

			for (int i = 0; i < Math.min(Farrayalpha.length, numIt + 1); i++) {
				for (int j = i; j < Math.min(Farrayalpha.length, numIt + 1);
					 j++) {
					mat.set(i, j, B.get(i, j));
					mat.set(j, i, B.get(i, j));

				}
			}

			double[] a = new double[mat.numRows()];
			Arrays.fill(a, 1);
			mat.setColumn(mat.numCols() - 1, 0, a);
			mat.setRow(mat.numRows() - 1, 0, a);
			mat.set(mat.numRows() - 1, mat.numCols() - 1, 0);

			SimpleMatrix rhs = new SimpleMatrix(mat.numRows(), 1);

			rhs.set(mat.numRows() - 1, 0, 1);

			try {
				SimpleMatrix DIIS = mat.solve(rhs);

				SimpleMatrix Fa = new SimpleMatrix(alphaDensity.numRows(), alphaDensity.numCols());
				SimpleMatrix Fb = new SimpleMatrix(betaDensity.numRows(), betaDensity.numCols());


				SimpleMatrix Da = new SimpleMatrix(alphaDensity.numRows(), alphaDensity.numCols());
				SimpleMatrix Db = new SimpleMatrix(betaDensity.numRows(), betaDensity.numCols());


				for (int i = 0; i < DIIS.getNumElements() - 1; i++) {
					Fa.plusi(Farrayalpha[i].scale(DIIS.get(i)));
					Da.plusi(Darrayalpha[i].scale(DIIS.get(i)));

					Fb.plusi(Farraybeta[i].scale(DIIS.get(i)));
					Db.plusi(Darraybeta[i].scale(DIIS.get(i)));
				}


				SimpleMatrix[] matrices1 = Utils.symEigen(Fa);

				SimpleMatrix[] matrices2 = Utils.symEigen(Fb);

				Ea = matrices1[1].diag();

				Eb = matrices2[1].diag();

				ca = matrices1[0].transpose();
				cb = matrices2[0].transpose();


				if (ca.get(0, 0) != ca.get(0, 0)) {

					matrices1 = Utils.symEigen(this.Fa);

					matrices2 = Utils.symEigen(this.Fb);

					Ea = matrices1[1].diag();

					Eb = matrices2[1].diag();

					ca = matrices1[0].transpose();
					cb = matrices2[0].transpose();


				}

				alphaDensity = calculateDensityMatrix(ca, nalpha);

				betaDensity = calculateDensityMatrix(cb, nbeta);


			} catch (Exception e) {
				SimpleMatrix[] matrices1 = Utils.symEigen(Fa);

				SimpleMatrix[] matrices2 = Utils.symEigen(Fb);

				Ea = matrices1[1].diag();

				Eb = matrices2[1].diag();

				ca = matrices1[0].transpose();
				cb = matrices2[0].transpose();

				alphaDensity = calculateDensityMatrix(ca, nalpha).scale(1 - damp)
						.plus(oldalphadensity.scale(damp));

				betaDensity = calculateDensityMatrix(cb, nbeta).scale(1 - damp)
						.plus(oldbetadensity.scale(damp));
			}

			if (numIt > 10000) {
				getRm().getLogger().error("unstable");
				System.exit(0);
			}

			numIt++;


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

		this.homo = Ea.get(nalpha - 1, 0);
		this.lumo = 0.001 * Math.round(Eb.get(nbeta, 0) * 1000);
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

			if (orbsOfAtom[j].length > 1) {//exclude hydrogen
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

		//System.err.println ("numit: " + numIt);
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
									nom.OneCenterERI(orbitals[i],
											orbitals[j],
											orbitals[k], orbitals[l]);
							e -= 0.5 * (alphaDensity.get(i, j) *
									alphaDensity.get(k, l) +
									betaDensity.get(i, j) *
											betaDensity.get(k, l)) *
									nom.OneCenterERI(orbitals[i],
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
							nom.beta(orbitals[i], orbitals[k]);
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
									* nom.getG(orbitals[i], orbitals[j],
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
	public SimpleMatrix alphaDensity() {
		return this.alphaDensity;
	}

	@Override
	public SimpleMatrix betaDensity() {
		return this.betaDensity;
	}

	@Override
	public Solution withNewAtoms(NDDOAtom[] newAtoms) {
		return new SolutionU(rm, newAtoms).compute();
	}

	@Override
	public SimpleMatrix densityMatrix() {
		return this.alphaDensity.plus(this.betaDensity);
	}
}
