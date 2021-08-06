package nddoparam;

import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import org.jblas.exceptions.LapackException;
import runcycle.input.RawMolecule;
import scf.Utils;

import java.util.ArrayList;

public class SolutionR extends Solution {
	private static final int[][][] TBRS =
			new int[][][]{new int[][]{new int[]{}},
					new int[][]{new int[]{0}},
					new int[][]{new int[]{1}, new int[]{0, 1}},
					new int[][]{new int[]{2}, new int[]{0, 2},
							new int[]{1, 2}, new int[]{0, 1, 2}},
					new int[][]{new int[]{3}, new int[]{0, 3}, new int[]{1, 3},
							new int[]{2, 3},
							new int[]{0, 1, 3}, new int[]{0, 2, 3},
							new int[]{1, 2, 3}, new int[]{0, 1, 2, 3}},
					new int[][]{new int[]{4}, new int[]{0, 4}, new int[]{1, 4},
							new int[]{2, 4}, new int[]{3, 4},
							new int[]{0, 1, 4}, new int[]{0, 2, 4},
							new int[]{0, 3, 4}, new int[]{1, 2, 4},
							new int[]{1, 3, 4}, new int[]{2, 3, 4},
							new int[]{0, 1, 2, 4}, new int[]{0, 1, 3, 4},
							new int[]{0, 2, 3, 4}, new int[]{1, 2, 3, 4},
							new int[]{0, 1, 2, 3, 4}},
					new int[][]{new int[]{5}, new int[]{0, 5}, new int[]{1, 5},
							new int[]{2, 5}, new int[]{3, 5}, new int[]{4, 5},
							new int[]{0, 1, 5}, new int[]{0, 2, 5},
							new int[]{0, 3, 5}, new int[]{0, 4, 5},
							new int[]{1, 2, 5}, new int[]{1, 3, 5},
							new int[]{1, 4, 5}, new int[]{2, 3, 5},
							new int[]{2, 4, 5}, new int[]{3, 4, 5},
							new int[]{0, 1, 2, 5}, new int[]{0, 1, 3, 5},
							new int[]{0, 1, 4, 5}, new int[]{0, 2, 3, 5},
							new int[]{0, 2, 4, 5}, new int[]{0, 3, 4, 5},
							new int[]{1, 2, 3, 5}, new int[]{1, 2, 4, 5},
							new int[]{1, 3, 4, 5}, new int[]{2, 3, 4, 5},
							new int[]{0, 1, 2, 3, 5}, new int[]{0, 1, 2, 4, 5},
							new int[]{0, 1, 3, 4, 5}, new int[]{0, 2, 3, 4, 5},
							new int[]{1, 2, 3, 4, 5},
							new int[]{0, 1, 2, 3, 4, 5}},
					new int[][]{new int[]{6}, new int[]{0, 6}, new int[]{1, 6},
							new int[]{2, 6}, new int[]{3, 6}, new int[]{4, 6},
							new int[]{5, 6}, new int[]{0, 1, 6},
							new int[]{0, 2, 6}, new int[]{0, 3, 6},
							new int[]{0, 4, 6}, new int[]{0, 5, 6},
							new int[]{1, 2, 6}, new int[]{1, 3, 6},
							new int[]{1, 4, 6}, new int[]{1, 5, 6},
							new int[]{2, 3, 6}, new int[]{2, 4, 6},
							new int[]{2, 5, 6}, new int[]{3, 4, 6},
							new int[]{3, 5, 6}, new int[]{4, 5, 6},
							new int[]{0, 1, 2, 6}, new int[]{0, 1, 3, 6},
							new int[]{0, 1, 4, 6}, new int[]{0, 1, 5, 6},
							new int[]{0, 2, 3, 6}, new int[]{0, 2, 4, 6},
							new int[]{0, 2, 5, 6}, new int[]{0, 3, 4, 6},
							new int[]{0, 3, 5, 6}, new int[]{0, 4, 5, 6},
							new int[]{1, 2, 3, 6}, new int[]{1, 2, 4, 6},
							new int[]{1, 2, 5, 6}, new int[]{1, 3, 4, 6},
							new int[]{1, 3, 5, 6}, new int[]{1, 4, 5, 6},
							new int[]{2, 3, 4, 6}, new int[]{2, 3, 5, 6},
							new int[]{2, 4, 5, 6}, new int[]{3, 4, 5, 6},
							new int[]{0, 1, 2, 3, 6}, new int[]{0, 1, 2, 4, 6},
							new int[]{0, 1, 2, 5, 6}, new int[]{0, 1, 3, 4, 6},
							new int[]{0, 1, 3, 5, 6}, new int[]{0, 1, 4, 5, 6},
							new int[]{0, 2, 3, 4, 6}, new int[]{0, 2, 3, 5, 6},
							new int[]{0, 2, 4, 5, 6}, new int[]{0, 3, 4, 5, 6},
							new int[]{1, 2, 3, 4, 6}, new int[]{1, 2, 3, 5, 6},
							new int[]{1, 2, 4, 5, 6}, new int[]{1, 3, 4, 5, 6},
							new int[]{2, 3, 4, 5, 6},
							new int[]{0, 1, 2, 3, 4, 6},
							new int[]{0, 1, 2, 3, 5, 6},
							new int[]{0, 1, 2, 4, 5, 6},
							new int[]{0, 1, 3, 4, 5, 6},
							new int[]{0, 2, 3, 4, 5, 6},
							new int[]{1, 2, 3, 4, 5, 6},
							new int[]{0, 1, 2, 3, 4, 5, 6}}};
	int integralEvaled, integralCached;
	public double[] integralArray;
	// H - core matrix, G = 2-electron matrix, F = fock matrix, C = coeffecient
	// matrix (transposed for easier reading), E = eigenvalues
	public DoubleMatrix C, F, G, E;
	private DoubleMatrix densityMatrix, B;
	private double[] Earray;
	private double[][][][] integralCache;
	private boolean[][][][] hasIntegralCache;

	public SolutionR(NDDOAtom[] atoms, RawMolecule rm) {
		super(atoms, rm);
		integralCache = new double[nOrbitals][nOrbitals]
				[nOrbitals][nOrbitals];
		hasIntegralCache = new boolean[nOrbitals][nOrbitals]
				[nOrbitals][nOrbitals];
	}

	private static DoubleMatrix commutator(DoubleMatrix F, DoubleMatrix D) {
		return F.mmul(D).sub(D.mmul(F));
	}

	private static DoubleMatrix removeElementsSquare(DoubleMatrix original,
													 int[] indices) {
		DoubleMatrix newarray = DoubleMatrix
				.zeros(original.rows - indices.length,
						original.rows - indices.length);

		ArrayList<Integer> array = new ArrayList<>();
		for (int i = 0; i < original.rows; i++) {
			array.add(i);
		}

		for (int i : indices) {
			array.remove(Integer.valueOf(i));
		}

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
													 int[] indices) {//get rid
		// of the rows given in indices and return downsized vector

		DoubleMatrix newarray =
				DoubleMatrix.zeros(original.rows - indices.length, 1);

		ArrayList<Integer> array = new ArrayList<>();
		for (int i = 0; i < original.rows; i++) {
			array.add(i);
		}

		for (int i : indices) {
			array.remove(Integer.valueOf(i));
		}

		int count = 0;

		for (int i : array) {
			newarray.put(count, original.get(i));
			count++;
		}

		return newarray;
	}

	private static DoubleMatrix addRows(DoubleMatrix original,
										int[] indices) { // add zero row at
		// indices

		DoubleMatrix newarray =
				DoubleMatrix.zeros(original.rows + indices.length, 1);

		ArrayList<Double> array = new ArrayList<>();

		for (double i : original.toArray()) {
			array.add(i);
		}

		for (int index : indices) {
			array.add(index, 0.0);
		}

		for (int i = 0; i < array.size(); i++) {
			newarray.put(i, array.get(i));
		}

		return newarray;
	}

	@Override
	public SolutionR compute() {
		StopWatch sw = new StopWatch();
		sw.start();

		integralArray = new double[rm.nIntegrals];

		int integralCount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							integralArray[integralCount] =
									aabb(j, l) - 0.5 * abab(j, l);
							integralCount++;
						}
					}

					for (int l : missingOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (m > -1) {
									if (atomOfOrb[l] == atomOfOrb[m]) {
										integralArray[integralCount] =
												aabc(j, l, m);
										integralCount++;
									}
								}

							}
						}
					}
				}
				else if (atomOfOrb[j] == atomOfOrb[k]) {
					double jkjk = abab(j, k);
					double jjkk = aabb(j, k);
					integralArray[integralCount] = 1.5 * jkjk - 0.5 * jjkk;
					integralCount++;

					for (int l : missingOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (m > -1) {
									if (atomOfOrb[l] == atomOfOrb[m]) {
										integralArray[integralCount] =
												abcd(j, k, l, m);
										integralCount++;
									}
								}
							}
						}
					}
				}
				else {
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
								if (m > -1) {
									integralArray[integralCount] =
											-0.5 * abcd(j, l, k, m);
									integralCount++;
								}
							}
						}
					}
				}
			}
		}

		System.out.println("integralEvals = " + integralEvaled);
		System.out.println("integralCached = " + integralCached);

		DoubleMatrix[] matrices = Utils.symEigen(H);
		E = matrices[1].diag();
		C = matrices[0].transpose();
		G = DoubleMatrix.zeros(C.rows, C.columns);
		F = H.dup();
		densityMatrix = calculateDensityMatrix(C);
		DoubleMatrix olddensity;

		DoubleMatrix[] Farray = new DoubleMatrix[8];
		DoubleMatrix[] Darray = new DoubleMatrix[8];
		Earray = new double[8];

		DoubleMatrix Bforediis = DoubleMatrix.zeros(8, 8);
		B = DoubleMatrix.zeros(8, 8);

		DoubleMatrix[] commutatorarray = new DoubleMatrix[8];

		int numIt = 0;
		double DIISError = 10;

		while (DIISError > 1E-11) {
			olddensity = densityMatrix.dup();
			integralCount = 0;

			// this entire block of code fills up the G matrix, and it calls
			// the integralarray to save time.

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {

						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								val += densityMatrix.get(l, l) *
										integralArray[integralCount];
								integralCount++;
							}
						}

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : missingOfAtom[atomOfOrb[j]]) {
									if (m > -1) {
										if (atomOfOrb[l] == atomOfOrb[m]) {
											val += densityMatrix.get(l, m) *
													integralArray[integralCount];
											integralCount++;
										}
									}

								}
							}
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						val += densityMatrix.get(j, k) *
								integralArray[integralCount];
						integralCount++;

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : missingOfAtom[atomOfOrb[j]]) {
									if (m > -1) {
										if (atomOfOrb[l] == atomOfOrb[m]) {
											val += densityMatrix.get(l, m) *
													integralArray[integralCount];
											integralCount++;
										}
									}

								}
							}
						}
					}
					else {
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : orbsOfAtom[atomOfOrb[k]]) {
									if (m > -1) {
										val += densityMatrix.get(l, m) *
												integralArray[integralCount];
										integralCount++;
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

			int ediisSize = Math.min(Farray.length + 1, numIt + 2);
			if (commutatorarray[Math.min(Farray.length - 1, numIt)].max() >
					0.01) {
				// if true do EDIIS else DIIS
				DoubleMatrix mat = DoubleMatrix.zeros(ediisSize, ediisSize);

				for (int i = 0; i < ediisSize - 1; i++) {
					for (int j = i; j < ediisSize - 1; j++) {
						mat.put(i, j, Bforediis.get(i, j));
						mat.put(j, i, Bforediis.get(i, j));
					}
				}

				mat.putColumn(mat.columns - 1, DoubleMatrix.ones(mat.rows, 1));
				mat.putRow(mat.rows - 1, DoubleMatrix.ones(mat.columns, 1));
				mat.put(mat.rows - 1, mat.columns - 1, 0);

				DoubleMatrix rhs = DoubleMatrix.ones(mat.rows, 1);
				for (int i = 0; i < ediisSize - 1; i++) {
					rhs.put(i, Earray[i]);
				}

				double bestE = 0;
				DoubleMatrix bestDIIS = null;
				try {
					int n = mat.rows - 2;
					for (int i = 0; i <= n; i++) {
						for (int[] tbr : TBRS[i]) {
							DoubleMatrix newmat =
									removeElementsSquare(mat, tbr);
							DoubleMatrix newrhs =
									removeElementsLinear(rhs, tbr);
							DoubleMatrix tempEdiis =
									addRows(Utils.solve(newmat, newrhs), tbr);
							tempEdiis.put(tempEdiis.rows - 1, 0);
							boolean nonNegative = !(tempEdiis.min() < 0);

							if (nonNegative) {

								double e = finde(tempEdiis);
								if (e < bestE) {
									bestE = e;
									bestDIIS = tempEdiis;
								}
							}
						}
					}
				} catch (LapackException ignored) {
				}

				DoubleMatrix finalDIIS = bestDIIS;
				DoubleMatrix F = DoubleMatrix.zeros(densityMatrix.rows,
						densityMatrix.columns);

				for (int i = 0; i < finalDIIS.length - 1; i++) {
					F = F.add(Farray[i].mmul(finalDIIS.get(i)));
				}


				matrices = Utils.symEigen(F);

				E = matrices[1].diag();

				C = matrices[0].transpose();

				if (C.get(0, 0) != C.get(0, 0)) {
					//System.err.println(
//							"NaN occurred at very much not DIIS iteration " +
//							numIt);
					//System.err.println("Exiting very much not DIIS for 1
					// Iteration...");

					matrices = Utils.symEigen(this.F);

					E = matrices[1].diag();

					C = matrices[0].transpose();


				}

				densityMatrix = calculateDensityMatrix(C);

			}
			else {
				DoubleMatrix mat = DoubleMatrix
						.zeros(ediisSize, ediisSize);

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
					DoubleMatrix DIIS = Utils.solve(mat, rhs);

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


					matrices = Utils.symEigen(F);

					E = matrices[1].diag();

					C = matrices[0].transpose();

					if (C.get(0, 0) != C.get(0, 0)) {
						//System.err.println(
//								"NaN occurred at DIIS iteration " + numIt);
						//
						//System.err.println("Exiting DIIS for 1 Iteration..
						// .");

						matrices = Utils.symEigen(this.F);

						E = matrices[1].diag();

						C = matrices[0].transpose();


					}

					densityMatrix = calculateDensityMatrix(C);
				} catch (Exception e) {
					matrices = Utils.symEigen(F);

					E = matrices[1].diag();

					C = matrices[0].transpose();

					double damp = 0.8;
					densityMatrix = calculateDensityMatrix(C).mmul(1 - damp)
							.add(olddensity.mmul(damp));
				}
			}

			numIt++;
		}

		findHF();
		findHomo();
		findDipole();

		System.out.println("sw.getTime()iamstupid = " + sw.getTime());
		return this;
	}

	private double aabb(int a, int b) {
		double aabb;
		if (hasIntegralCache[a][a][b][b]) {
			aabb = integralCache[a][a][b][b];
			integralCached++;
		}
		else {
			aabb = NDDO6G.OneCenterERI(
					orbitals[a],
					orbitals[a],
					orbitals[b],
					orbitals[b]);
			integralCache[a][a][b][b] = aabb;
			integralCache[b][b][a][a] = aabb;

			hasIntegralCache[a][a][b][b] = true;
			hasIntegralCache[b][b][a][a] = true;
			integralEvaled++;
		}
		return aabb;
	}

	private double aabc(int a, int b, int c) {
		double aabc;
		if (hasIntegralCache[a][a][b][c]) {
			aabc = integralCache[a][a][b][c];
			integralCached++;
		}
		else {
			aabc = NDDO6G.getG(orbitals[a],
					orbitals[a],
					orbitals[b],
					orbitals[c]);
			integralCache[a][a][b][c] = aabc;
			integralCache[a][a][c][b] = aabc;
			integralCache[b][c][a][a] = aabc;
			integralCache[c][b][a][a] = aabc;
			hasIntegralCache[a][a][b][c] =
					true;
			hasIntegralCache[a][a][c][b] =
					true;
			hasIntegralCache[b][c][a][a] =
					true;
			hasIntegralCache[c][b][a][a] =
					true;
			integralEvaled++;
		}
		return aabc;
	}

	private double abab(int a, int b) {
		double abab;
		if (hasIntegralCache[a][b][a][b]) {
			abab = integralCache[a][b][a][b];
			integralCached++;
		}
		else {
			abab = NDDO6G.OneCenterERI(
					orbitals[a],
					orbitals[b],
					orbitals[a],
					orbitals[b]);
			integralCache[a][b][a][b] = abab;
			integralCache[a][b][b][a] = abab;
			integralCache[b][a][a][b] = abab;
			integralCache[b][a][b][a] = abab;

			hasIntegralCache[a][b][a][b] = true;
			hasIntegralCache[a][b][b][a] = true;
			hasIntegralCache[b][a][a][b] = true;
			hasIntegralCache[b][a][b][a] = true;
			integralEvaled++;
		}
		return abab;
	}

	private double abcd(int a, int b, int c, int d) {
		double abcd;
		if (hasIntegralCache[a][b][c][d]) {
			abcd = integralCache[a][b][c][d];
			integralCached++;
		}
		else {
			abcd = NDDO6G.OneCenterERI(
					orbitals[a],
					orbitals[b],
					orbitals[c],
					orbitals[d]);
			integralCache[a][b][c][d] = abcd;
			integralCache[a][b][d][c] = abcd;
			integralCache[b][a][c][d] = abcd;
			integralCache[b][a][d][c] = abcd;
			integralCache[c][d][a][b] = abcd;
			integralCache[c][d][b][a] = abcd;
			integralCache[d][c][a][b] = abcd;
			integralCache[d][c][b][a] = abcd;

			hasIntegralCache[a][b][c][d] = true;
			hasIntegralCache[a][b][d][c] = true;
			hasIntegralCache[b][a][c][d] = true;
			hasIntegralCache[b][a][d][c] = true;
			hasIntegralCache[c][d][a][b] = true;
			hasIntegralCache[c][d][b][a] = true;
			hasIntegralCache[d][c][a][b] = true;
			hasIntegralCache[d][c][b][a] = true;
			integralEvaled++;
		}
		return abcd;
	}

	@SuppressWarnings("DuplicatedCode")
	protected int findNIntegrals() {
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

	private void findHF() {
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				energy += 0.5 * densityMatrix.get(j, k) *
						(H.get(j, k) + F.get(j, k));
			}
		}

		double heat = 0;
		for (int j = 0; j < atoms.length; j++) {
			heat += atoms[j].getHeat() - atoms[j].getEisol();
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

			if (orbsOfAtom[j][1] != -1) {//exclude hydrogen
				hybridip[0] = hybridip[0] - 2.5416 * 2 * atoms[j].D1 *
						densityMatrix.get(orbsOfAtom[j][0],
								orbsOfAtom[j][1]);
				hybridip[1] = hybridip[1] - 2.5416 * 2 * atoms[j].D1 *
						densityMatrix.get(orbsOfAtom[j][0],
								orbsOfAtom[j][2]);
				hybridip[2] = hybridip[2] - 2.5416 * 2 * atoms[j].D1 *
						densityMatrix.get(orbsOfAtom[j][0],
								orbsOfAtom[j][3]);
			}
		}

		dipoletot = new double[]{chargedip[0] + hybridip[0],
				chargedip[1] + hybridip[1],
				chargedip[2] + hybridip[2]};

		dipole = Math.sqrt(
				dipoletot[0] * dipoletot[0] + dipoletot[1] * dipoletot[1] +
						dipoletot[2] * dipoletot[2]);
	}

	private void findHomo() {
		if (nElectrons > 0) homo = E.get(nElectrons / 2 - 1, 0);
		else homo = 0;
		lumo = E.get(nElectrons / 2, 0);
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
}
