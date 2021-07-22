package nddoparam;

import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.jblas.Solve;
import runcycle.input.RawMolecule;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


public class SolutionNew extends Solution {

	public double[] integralArray;
	public DoubleMatrix C, F, G, E;
	//H - core matrix, G = 2-electron matrix, F = fock matrix, C = coeffecient
	// matrix (transposed for easier reading), E = eigenvalues
	private DoubleMatrix densityMatrix, B;
	private double[] Earray;

	public SolutionNew(NDDOAtom[] atoms, int charge) {
		super(atoms, charge);

		StopWatch sw = new StopWatch();

		sw.start();

		int size = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {

					for (int l : orbitalIndices[atomNumber[j]]) {
						if (l > -1) {
							size++;
						}
					}

					for (int l : missingIndex[atomNumber[j]]) {
						if (l > -1) {
							for (int m : missingIndex[atomNumber[j]]) {
								if (m > -1) {
									if (atomNumber[l] == atomNumber[m]) {
										size++;
									}
								}

							}
						}
					}
				}
				else if (atomNumber[j] == atomNumber[k]) {
					size++;

					for (int l : missingIndex[atomNumber[j]]) {
						if (l > -1) {
							for (int m : missingIndex[atomNumber[j]]) {
								if (m > -1) {
									if (atomNumber[l] == atomNumber[m]) {
										size++;
									}
								}

							}
						}
					}
				}
				else {
					for (int l : orbitalIndices[atomNumber[j]]) {
						if (l > -1) {
							for (int m : orbitalIndices[atomNumber[k]]) {
								if (m > -1) {
									size++;
								}
							}
						}
					}
				}
			}
		}
		integralArray = new double[size];
		//The idea of the integralarray is to simply store all the integrals
		// in order they are called. It's basically my way of avoiding
		// having to perform a Yoshemine sort.
		// TODO re-implement HashMap
		int integralcount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) { // case 1
					for (int l : orbitalIndices[atomNumber[j]]) {
						if (l > -1) {
							integralArray[integralcount] =
									(NDDO6G.OneCenterERI(orbitals[j],
											orbitals[j],
											orbitals[l], orbitals[l]) - 0.5 *
											NDDO6G.OneCenterERI(orbitals[j],
													orbitals[l],
													orbitals[j], orbitals[l]));
							integralcount++;
						}
					}

					for (int l : missingIndex[atomNumber[j]]) {
						if (l > -1) {
							for (int m : missingIndex[atomNumber[j]]) {
								if (m > -1) {
									if (atomNumber[l] == atomNumber[m]) {
										integralArray[integralcount] =
												(NDDO6G.getG(orbitals[j],
														orbitals[j],
														orbitals[l],
														orbitals[m]));
										integralcount++;
									}
								}

							}
						}
					}
				}
				else if (atomNumber[j] == atomNumber[k]) { // case 2
					integralArray[integralcount] = (1.5 *
							NDDO6G.OneCenterERI(orbitals[j], orbitals[k],
									orbitals[j],
									orbitals[k]) - 0.5 *
							NDDO6G.OneCenterERI(orbitals[j], orbitals[j],
									orbitals[k],
									orbitals[k]));
					integralcount++;
					for (int l : missingIndex[atomNumber[j]]) {
						if (l > -1) {
							for (int m : missingIndex[atomNumber[j]]) {
								if (m > -1) {
									if (atomNumber[l] == atomNumber[m]) {
										integralArray[integralcount] =
												(NDDO6G.getG(orbitals[j],
														orbitals[k],
														orbitals[l],
														orbitals[m]));
										integralcount++;
									}
								}
							}
						}
					}
				}
				else { // case 3
					for (int l : orbitalIndices[atomNumber[j]]) {
						if (l > -1) {
							for (int m : orbitalIndices[atomNumber[k]]) {
								if (m > -1) {
									integralArray[integralcount] = (-0.5 *
											NDDO6G.getG(orbitals[j],
													orbitals[l],
													orbitals[k], orbitals[m]));
									integralcount++;
								}
							}
						}
					}
				}
			}
		}

		DoubleMatrix[] matrices = Eigen.symmetricEigenvectors(H);
		E = matrices[1].diag();

		C = matrices[0].transpose();

		G = DoubleMatrix.zeros(C.rows, C.columns);

		densityMatrix = calculateDensityMatrix(C);

		DoubleMatrix olddensity = DoubleMatrix.zeros(C.rows, C.columns);

		F = H.dup();

		DoubleMatrix[] Farray = new DoubleMatrix[8];
		DoubleMatrix[] Darray = new DoubleMatrix[8];
		Earray = new double[8];

		DoubleMatrix Bforediis = DoubleMatrix.zeros(8, 8);

		B = DoubleMatrix.zeros(8, 8);

		DoubleMatrix[] commutatorarray = new DoubleMatrix[8];


		int numIt = 0;

		double DIISError = 10;
		while (DIISError > 1E-10) {
			System.out.println("numIt = " + numIt);
			olddensity = densityMatrix.dup();

			integralcount = 0;

			//this entire block of code fills up the G matrix, and it calls the
			// integralarray to save time.

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {

						for (int l : orbitalIndices[atomNumber[j]]) {
							if (l > -1) {
								val += densityMatrix.get(l, l) *
										integralArray[integralcount];
								integralcount++;
							}
						}

						for (int l : missingIndex[atomNumber[j]]) {
							if (l > -1) {
								for (int m : missingIndex[atomNumber[j]]) {
									if (m > -1) {
										if (atomNumber[l] == atomNumber[m]) {
											val += densityMatrix.get(l, m) *
													integralArray[integralcount];
											integralcount++;
										}
									}

								}
							}
						}
					}
					else if (atomNumber[j] == atomNumber[k]) {
						val += densityMatrix.get(j, k) *
								integralArray[integralcount];
						integralcount++;

						for (int l : missingIndex[atomNumber[j]]) {
							if (l > -1) {
								for (int m : missingIndex[atomNumber[j]]) {
									if (m > -1) {
										if (atomNumber[l] == atomNumber[m]) {
											val += densityMatrix.get(l, m) *
													integralArray[integralcount];
											integralcount++;
										}
									}

								}
							}
						}
					}
					else {
						for (int l : orbitalIndices[atomNumber[j]]) {
							if (l > -1) {
								for (int m : orbitalIndices[atomNumber[k]]) {
									if (m > -1) {
										val += densityMatrix.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}
							}
						}
					}

					G.put(j, k, val);
					G.put(k, j, val);
				}
			}

			F = H.dup().add(G);

			if (numIt < Farray.length) {

				Farray[numIt] = F.dup();

				Darray[numIt] = densityMatrix.dup();
				Earray[numIt] = -0.5 * (H.mmul(densityMatrix)).diag().sum();
				commutatorarray[numIt] =
						commutator(F.dup(), densityMatrix.dup());
				DIISError = commutatorarray[numIt].norm2();

				for (int i = 0; i <= numIt; i++) {

					double product =
							(commutatorarray[numIt]
									.mmul(commutatorarray[i].transpose()))
									.diag().sum();
					B.put(i, numIt, product);
					B.put(numIt, i, product);

					product = 0.5 *
							((Farray[i].mmul(Darray[numIt])).diag().sum() +
									(Farray[numIt].mmul(Darray[i])).diag()
											.sum());

					Bforediis.put(i, numIt, product);
					Bforediis.put(numIt, i, product);
				}
			}
			else {

				for (int i = 0; i < Farray.length - 1; i++) {

					Farray[i] = Farray[i + 1].dup();
					Darray[i] = Darray[i + 1].dup();
					commutatorarray[i] = commutatorarray[i + 1].dup();
					Earray[i] = Earray[i + 1];
				}

				Farray[Farray.length - 1] = F.dup();
				Darray[Darray.length - 1] = densityMatrix.dup();
				Earray[Darray.length - 1] =
						-0.5 * (H.mmul(densityMatrix)).diag().sum();
				commutatorarray[Darray.length - 1] =
						commutator(F.dup(), densityMatrix.dup());
				DIISError = commutatorarray[Darray.length - 1].norm2();

				// B is dy/dx sort of, make dy/dx 0
				DoubleMatrix newB = DoubleMatrix.zeros(8, 8);

				DoubleMatrix newBforediis = DoubleMatrix.zeros(8, 8);

				for (int i = 0; i < Farray.length - 1; i++) {
					for (int j = i; j < Farray.length - 1; j++) {
						newB.put(i, j, B.get(i + 1, j + 1));
						newB.put(j, i, B.get(i + 1, j + 1));
						newBforediis.put(i, j, Bforediis.get(i + 1, j + 1));
						newBforediis.put(j, i, Bforediis.get(i + 1, j + 1));
					}
				}

				for (int i = 0; i < Farray.length; i++) {

					double product =
							commutatorarray[Farray.length - 1].transpose()
									.mmul(commutatorarray[i]).diag().sum();
					newB.put(i, Farray.length - 1, product);
					newB.put(Farray.length - 1, i, product);

					product = 0.5 *
							((Farray[i].mmul(Darray[Farray.length - 1])).diag()
									.sum() +
									(Farray[Farray.length - 1].mmul(Darray[i]))
											.diag()
											.sum());
					newBforediis.put(i, Farray.length - 1, product);
					newBforediis.put(Farray.length - 1, i, product);
				}

				B = newB.dup();

				Bforediis = newBforediis.dup();
			}


			//System.err.println ("DIIS Error: " + DIISError);


			if (commutatorarray[Math.min(Farray.length - 1, numIt)].max() >
					0.01) {
				// if true do EDIIS else DIIS

				DoubleMatrix mat = DoubleMatrix
						.zeros(Math.min(Farray.length + 1, numIt + 2),
								Math.min(Farray.length + 1, numIt + 2));

				for (int i = 0; i < Math.min(Farray.length, numIt + 1); i++) {
					for (int j = i; j < Math.min(Farray.length, numIt + 1);
						 j++) {
						mat.put(i, j, Bforediis.get(i, j));
						mat.put(j, i, Bforediis.get(i, j));

					}
				}


				mat.putColumn(mat.columns - 1, DoubleMatrix.ones(mat.rows, 1));

				mat.putRow(mat.rows - 1, DoubleMatrix.ones(mat.columns, 1));

				mat.put(mat.rows - 1, mat.columns - 1, 0);

				DoubleMatrix rhs = DoubleMatrix.ones(mat.rows, 1);
				System.out.println("Earray = " + Arrays.toString(Earray));
				for (int i = 0; i < Math.min(Farray.length, numIt + 1); i++) {
					rhs.put(i, Earray[i]);
				}
//
//
//				EdiisTry firstTry = null;
//				DoubleMatrix DIIS = Solve.solve(mat, rhs);
//
//				DIIS = DIIS.put(DIIS.rows - 1, 0);
//
//				boolean nonNegative = !(DIIS.min() < 0);
//
//				if (nonNegative) {
//					double e = finde(DIIS);
//					firstTry = new EdiisTry(DIIS, e);
//				}
				EdiisTry bestDIIS =
						findBestEdiis(mat, rhs, new ArrayList<>(8), null);
				System.out.println("bestDIIS = " + bestDIIS.attempt);

				DoubleMatrix finalDIIS = bestDIIS.attempt;

				DoubleMatrix F = DoubleMatrix.zeros(densityMatrix.rows,
						densityMatrix.columns);

				for (int i = 0; i < finalDIIS.length - 1; i++) {
					F = F.add(Farray[i].mmul(finalDIIS.get(i)));
				}


				matrices = Eigen.symmetricEigenvectors(F);

				E = matrices[1].diag();

				C = matrices[0].transpose();

				if (C.get(0, 0) != C.get(0, 0)) {
					//System.err.println(
//							"NaN occurred at very much not DIIS iteration " +
//							numIt);
					//System.err.println("Exiting very much not DIIS for 1
					// Iteration...");

					matrices = Eigen.symmetricEigenvectors(this.F);

					E = matrices[1].diag();

					C = matrices[0].transpose();


				}

				densityMatrix = calculateDensityMatrix(C);

				System.out.println();
			}
			else {

				DoubleMatrix mat = DoubleMatrix
						.zeros(Math.min(Farray.length + 1, numIt + 2),
								Math.min(Farray.length + 1, numIt + 2));

				for (int i = 0; i < Math.min(Farray.length, numIt + 1); i++) {
					for (int j = i; j < Math.min(Farray.length, numIt + 1);
						 j++) {
						mat.put(i, j, B.get(i, j));
						mat.put(j, i, B.get(i, j));

					}
				}


				mat.putColumn(mat.columns - 1, DoubleMatrix.ones(mat.rows, 1));

				mat.putRow(mat.rows - 1, DoubleMatrix.ones(mat.columns, 1));

				mat.put(mat.rows - 1, mat.columns - 1, 0);

				DoubleMatrix rhs = DoubleMatrix.zeros(mat.rows, 1);

				rhs.put(mat.rows - 1, 0, 1);

				try {
					DoubleMatrix DIIS = Solve.solve(mat, rhs);

					DoubleMatrix F =
							DoubleMatrix.zeros(densityMatrix.rows,
									densityMatrix.columns);

					DoubleMatrix D =
							DoubleMatrix.zeros(densityMatrix.rows,
									densityMatrix.columns);


					for (int i = 0; i < DIIS.length - 1; i++) {
						F = F.add(Farray[i].mmul(DIIS.get(i)));
						D = D.add(Darray[i].mmul(DIIS.get(i)));
					}


					matrices = Eigen.symmetricEigenvectors(F);

					E = matrices[1].diag();

					C = matrices[0].transpose();

					if (C.get(0, 0) != C.get(0, 0)) {
						//System.err.println(
//								"NaN occurred at DIIS iteration " + numIt);
						//
						//System.err.println("Exiting DIIS for 1 Iteration..
						// .");

						matrices = Eigen.symmetricEigenvectors(this.F);

						E = matrices[1].diag();

						C = matrices[0].transpose();


					}

					densityMatrix = calculateDensityMatrix(C);
				} catch (Exception e) {
					matrices = Eigen.symmetricEigenvectors(F);

					E = matrices[1].diag();

					C = matrices[0].transpose();

					densityMatrix = calculateDensityMatrix(C).mmul(1 - damp)
							.add(olddensity.mmul(damp));
				}
			}


			numIt++;

		}

		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * densityMatrix.get(j, k) *
						(H.get(j, k) + F.get(j, k));
			}
		}

		double heat = 0;

		for (int j = 0; j < atoms.length; j++) {
			heat += atoms[j].getHeat() - atoms[j].getEisol();
			for (int k = j + 1; k < atoms.length; k++) {
				e += atoms[j].crf(atoms[k]);
			}
		}

		energy = e;

		if (energy != energy) System.err.println(densityMatrix.getRow(0));
		heat += e;

		this.hf = heat / 4.3363E-2;

		if (nElectrons > 0) {
			this.homo = E.get(nElectrons / 2 - 1, 0);
		}
		else {
			this.homo = 0;
		}

		this.lumo = E.get(nElectrons / 2, 0);

		double[] populations = new double[atoms.length];

		for (int j = 0; j < atoms.length; j++) {
			double sum = 0;
			for (int k : orbitalIndices[j]) {
				if (k > -1) {
					sum += densityMatrix.get(k, k);
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

			if (orbitalIndices[j][1] != -1) {//exclude hydrogen
				hybridip[0] = hybridip[0] - 2.5416 * 2 * atoms[j].D1 *
						densityMatrix.get(orbitalIndices[j][0],
								orbitalIndices[j][1]);
				hybridip[1] = hybridip[1] - 2.5416 * 2 * atoms[j].D1 *
						densityMatrix.get(orbitalIndices[j][0],
								orbitalIndices[j][2]);
				hybridip[2] = hybridip[2] - 2.5416 * 2 * atoms[j].D1 *
						densityMatrix.get(orbitalIndices[j][0],
								orbitalIndices[j][3]);
			}
		}


		dipoletot = new double[]{chargedip[0] + hybridip[0],
				chargedip[1] + hybridip[1],
				chargedip[2] + hybridip[2]};


		dipole = Math.sqrt(
				dipoletot[0] * dipoletot[0] + dipoletot[1] * dipoletot[1] +
						dipoletot[2] * dipoletot[2]);


	}

	private static DoubleMatrix commutator(DoubleMatrix F, DoubleMatrix D) {
		return F.mmul(D).sub(D.mmul(F));
	}

	private static DoubleMatrix removeElementsSquare(DoubleMatrix original,
													 List<Integer> indices,
													 List<Integer> array) {
		// remove rows and columns specified in indices
		// and return downsized square matrix

		int rows = original.rows - indices.size();
		DoubleMatrix newarray = DoubleMatrix.zeros(rows, rows);

		int count = 0;

		for (int i : array) {
			int count1 = 0;
			for (int j : array) {
				newarray.put(count, count1, original.get(i, j));
				count1++;
			}

			count++;
		}
		return newarray;
	}

	private static DoubleMatrix removeElementsLinear(DoubleMatrix original,
													 List<Integer> indices,
													 List<Integer> array) {
		//get rid
		// of the rows given in indices and return downsized vector

		DoubleMatrix newarray =
				DoubleMatrix.zeros(original.rows - indices.size(), 1);

		int count = 0;

		for (int i : array) {
			newarray.put(count, original.get(i));

			count++;
		}

		return newarray;
	}

	private static ArrayList<Integer> getComplement(DoubleMatrix original,
													List<Integer> indices) {
		ArrayList<Integer> array = new ArrayList<>();
		List<Integer> newIndices = new ArrayList<>(indices);
		for (int i = 0; i < original.rows; i++) {
			boolean in = true;
			for (int j = 0; j < newIndices.size(); j++) {
				if (newIndices.get(j) == i) {
					in = false;
					newIndices.remove(j);
					break;
				}
			}
			if (in) array.add(i);
		}
		return array;
	}

	private static DoubleMatrix addRows(DoubleMatrix original,
										List<Integer> indices) {
		// add zero row at indices
		DoubleMatrix res =
				DoubleMatrix.zeros(original.rows + indices.size(), 1);
		ArrayList<Double> padded = new ArrayList<>();
		double[] old = original.toArray();

		int oi = 0;
		for (int i = 0; i < 9; i++) {
			boolean skip = false;
			for (int index : indices) {
				if (index == i) {
					padded.add(0.0);
					skip = true;
					break;
				}
			}
			if (oi < old.length && !skip) {
				padded.add(old[oi]);
				oi++;
			}
		}

		for (int i = 0; i < padded.size(); i++) {
			res.put(i, padded.get(i));
		}

		return res;
	}

	private double finde(DoubleMatrix ediis) {
		double e = 0;

		for (int a = 0; a < ediis.length - 1; a++) {
			e -= Earray[a] * ediis.get(a);
		}

		for (int a = 0; a < ediis.length - 1; a++) {
			for (int b = 0; b < ediis.length - 1; b++) {
				e += 0.5 * ediis.get(a) * ediis.get(b) * 0.5 *
						B.get(a, b);
			}
		}
		return e;
	}

	private EdiisTry findBestEdiis(DoubleMatrix mat,
								   DoubleMatrix rhs,
								   List<Integer> tbrList,
								   EdiisTry bestEdiis) {
		// tbr stands for toBeRemoved
		int n = mat.rows - 1;
		if (bestEdiis != null && n-1==0) return bestEdiis;
		System.out.println("n = " + n);
		System.out.println("tbrList = " + tbrList);

		List<Integer> array = getComplement(mat, tbrList);
		DoubleMatrix smallmat = removeElementsSquare(mat.dup(), tbrList,
				array);
		DoubleMatrix smallrhs = removeElementsLinear(rhs.dup(), tbrList,
				array);
		System.out.println("smallmat = " + mat);
		System.out.println("smallrhs = " + rhs);
		DoubleMatrix attemptRaw = Solve.solve(smallmat, smallrhs);
		DoubleMatrix attempt = addRows(attemptRaw, tbrList);

		boolean nonNegative = !(attempt.min() < 0);

		// see if this try is better than best try
		if (nonNegative) {
			double e = finde(attempt);
			EdiisTry ediisTry = new EdiisTry(attempt, e);

			if (bestEdiis == null || ediisTry.e < bestEdiis.e) {
				bestEdiis = ediisTry;
			}
		}

		int size = tbrList.size();
		for (int i = size == 0 ? 0 : tbrList.get(size - 1) + 1; i < n; i++) {
			List<Integer> newTbrList = new ArrayList<>(tbrList);
			newTbrList.add(i);

			EdiisTry bestEdiisFurtherDown =
					findBestEdiis(mat, rhs, newTbrList, bestEdiis);

			// if a best ediis further down the tree is better than the
			// best one so far
			if (bestEdiis == null || bestEdiisFurtherDown.e < bestEdiis.e) {
				bestEdiis = bestEdiisFurtherDown;
			}
		}

		return bestEdiis;
	}

	@Override
	public SolutionNew setRm(RawMolecule rm) {
		this.rm = rm;
		return this;
	}

	private DoubleMatrix calculateDensityMatrix(
			DoubleMatrix c) {//density matrix construction by definition.
		DoubleMatrix densityMatrix = DoubleMatrix.zeros(orbitals.length,
				orbitals.length);
		for (int i = 0; i < orbitals.length; i++) {
			for (int j = 0; j < orbitals.length; j++) {
				double sum = 0;
				int count = nElectrons;
				int counter = -1;
				while (count > 0) {
					counter++;
					sum += 2 * c.get(counter, i) * c.get(counter, j);
					count -= 2;
				}
				densityMatrix.put(i, j, sum);
			}
		}
		return densityMatrix;
	}

	public DoubleMatrix getE() {
		return E;
	}

	@Override
	public DoubleMatrix alphaDensity() {
		return this.densityMatrix.mmul(0.5);
	}

	@Override
	public DoubleMatrix betaDensity() {
		return this.densityMatrix.mmul(0.5);
	}

	@Override
	public DoubleMatrix densityMatrix() {
		return this.densityMatrix;
	}

	class EdiisTry {
		DoubleMatrix attempt;
		double e;

		public EdiisTry(DoubleMatrix attempt, double e) {
			this.attempt = attempt;
			this.e = e;
		}
	}
}
