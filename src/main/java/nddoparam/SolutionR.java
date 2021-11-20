package nddoparam;

import org.ejml.data.SingularMatrixException;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.simple.SimpleMatrix;
import runcycle.input.RawMolecule;
import tools.Utils;

import java.util.ArrayList;
import java.util.Arrays;

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

	public double[] integralArray;
	public SimpleMatrix C, F, G, E;
	// H - core matrix, G = 2-electron matrix, F = fock matrix, C = coeffecient
	// matrix (transposed for easier reading), E = eigenvalues
	private SimpleMatrix densityMatrix, B;
	private double[] Earray;

	public SolutionR(NDDOAtom[] atoms, RawMolecule rm) {
		super(atoms, rm);
	}

	@Override
	public SolutionR compute() {
		integralArray = new double[rm.nIntegrals];

		int integralcount = 0;
		for (int j = 0; j < orbitals.length; j++) {
			for (int k = j; k < orbitals.length; k++) {
				if (j == k) { // case 1
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							integralArray[integralcount] =
									(NDDO6G.OneCenterERI(orbitals[j],
											orbitals[j],
											orbitals[l],
											orbitals[l]) - 0.5 *
											NDDO6G.OneCenterERI(orbitals[j],
													orbitals[l],
													orbitals[j],
													orbitals[l]));
							integralcount++;
						}
					}

					for (int l : missingOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (m > -1) {
									if (atomOfOrb[l] == atomOfOrb[m]) {
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
				else if (atomOfOrb[j] == atomOfOrb[k]) { // case 2
					integralArray[integralcount] = (1.5 *
							NDDO6G.OneCenterERI(orbitals[j], orbitals[k],
									orbitals[j],
									orbitals[k]) - 0.5 *
							NDDO6G.OneCenterERI(orbitals[j], orbitals[j],
									orbitals[k],
									orbitals[k]));
					integralcount++;
					for (int l : missingOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : missingOfAtom[atomOfOrb[j]]) {
								if (m > -1) {
									if (atomOfOrb[l] == atomOfOrb[m]) {
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
					for (int l : orbsOfAtom[atomOfOrb[j]]) {
						if (l > -1) {
							for (int m : orbsOfAtom[atomOfOrb[k]]) {
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

		SimpleMatrix[] matrices = Utils.symEigen(H);
		E = matrices[1].diag();
		C = matrices[0].transpose();
		G = new SimpleMatrix(C.numRows(), C.numCols());
		F = H.copy();
		densityMatrix = calculateDensityMatrix(C);
		SimpleMatrix olddensity;

		SimpleMatrix[] Farray = new SimpleMatrix[8];
		SimpleMatrix[] Darray = new SimpleMatrix[8];
		Earray = new double[8];

		SimpleMatrix Bforediis = new SimpleMatrix(8, 8);
		B = new SimpleMatrix(8, 8);

		SimpleMatrix[] commutatorarray = new SimpleMatrix[8];

		int numIt = 0;
		double DIISError = 10;

		while (DIISError > 1E-11) {
			olddensity = densityMatrix.copy();
			integralcount = 0;

			// this entire block of code fills up the G matrix, and it calls
			// the integralarray to save time.

			for (int j = 0; j < orbitals.length; j++) {
				for (int k = j; k < orbitals.length; k++) {
					double val = 0;
					if (j == k) {

						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								val += densityMatrix.get(l, l) *
										integralArray[integralcount];
								integralcount++;
							}
						}

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : missingOfAtom[atomOfOrb[j]]) {
									if (m > -1) {
										if (atomOfOrb[l] == atomOfOrb[m]) {
											val += densityMatrix.get(l, m) *
													integralArray[integralcount];
											integralcount++;
										}
									}

								}
							}
						}
					}
					else if (atomOfOrb[j] == atomOfOrb[k]) {
						val += densityMatrix.get(j, k) *
								integralArray[integralcount];
						integralcount++;

						for (int l : missingOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : missingOfAtom[atomOfOrb[j]]) {
									if (m > -1) {
										if (atomOfOrb[l] == atomOfOrb[m]) {
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
						for (int l : orbsOfAtom[atomOfOrb[j]]) {
							if (l > -1) {
								for (int m : orbsOfAtom[atomOfOrb[k]]) {
									if (m > -1) {
										val += densityMatrix.get(l, m) *
												integralArray[integralcount];
										integralcount++;
									}
								}
							}
						}
					}

					G.set(j, k, val);
					G.set(k, j, val);
				}
			}

			F = H.copy().plus(G);
			if (numIt < Farray.length) {

				Farray[numIt] = F.copy();

				Darray[numIt] = densityMatrix.copy();
				Earray[numIt] =
						-0.5 * (H.mult(densityMatrix)).diag().elementSum();
				commutatorarray[numIt] =
						commutator(F.copy(), densityMatrix.copy());
				DIISError = commutatorarray[numIt].normF();

				for (int i = 0; i <= numIt; i++) {

					double product =
							(commutatorarray[numIt]
									.mult(commutatorarray[i].transpose()))
									.diag().elementSum();
					B.set(i, numIt, product);
					B.set(numIt, i, product);

					product = 0.5 *
							((Farray[i].mult(Darray[numIt])).diag()
									.elementSum() +
									(Farray[numIt].mult(Darray[i])).diag()
											.elementSum());

					Bforediis.set(i, numIt, product);
					Bforediis.set(numIt, i, product);
				}
			}
			else {

				for (int i = 0; i < Farray.length - 1; i++) {

					Farray[i] = Farray[i + 1].copy();
					Darray[i] = Darray[i + 1].copy();
					commutatorarray[i] = commutatorarray[i + 1].copy();
					Earray[i] = Earray[i + 1];
				}

				Farray[Farray.length - 1] = F.copy();
				Darray[Darray.length - 1] = densityMatrix.copy();
				Earray[Darray.length - 1] =
						-0.5 * (H.mult(densityMatrix)).diag().elementSum();
				commutatorarray[Darray.length - 1] =
						commutator(F.copy(), densityMatrix.copy());
				DIISError = commutatorarray[Darray.length - 1].normF();

				// B is dy/dx sort of, make dy/dx 0
				SimpleMatrix newB = new SimpleMatrix(8, 8);

				SimpleMatrix newBforediis = new SimpleMatrix(8, 8);

				for (int i = 0; i < Farray.length - 1; i++) {
					for (int j = i; j < Farray.length - 1; j++) {
						newB.set(i, j, B.get(i + 1, j + 1));
						newB.set(j, i, B.get(i + 1, j + 1));
						newBforediis.set(i, j, Bforediis.get(i + 1, j + 1));
						newBforediis.set(j, i, Bforediis.get(i + 1, j + 1));
					}
				}

				for (int i = 0; i < Farray.length; i++) {

					double product =
							commutatorarray[Farray.length - 1].transpose()
									.mult(commutatorarray[i]).diag()
									.elementSum();
					newB.set(i, Farray.length - 1, product);
					newB.set(Farray.length - 1, i, product);

					product = 0.5 *
							((Farray[i].mult(Darray[Farray.length - 1])).diag()
									.elementSum() +
									(Farray[Farray.length - 1].mult(Darray[i]))
											.diag()
											.elementSum());
					newBforediis.set(i, Farray.length - 1, product);
					newBforediis.set(Farray.length - 1, i, product);
				}

				B = newB.copy();

				Bforediis = newBforediis.copy();
			}

			int ediisSize = Math.min(Farray.length + 1, numIt + 2);
			if (CommonOps_DDRM.elementMax(
					commutatorarray[Math.min(Farray.length - 1, numIt)]
							.getDDRM()) > 0.01) {
				// if true do EDIIS else DIIS
				SimpleMatrix mat = new SimpleMatrix(ediisSize, ediisSize);

				for (int i = 0; i < ediisSize - 1; i++) {
					for (int j = i; j < ediisSize - 1; j++) {
						mat.set(i, j, Bforediis.get(i, j));
						mat.set(j, i, Bforediis.get(i, j));
					}
				}

				double[] row = new double[mat.numRows()];
				double[] col = new double[mat.numCols()];
				Arrays.fill(row, 1);
				Arrays.fill(col, 1);
				mat.setColumn(mat.numCols() - 1, 0, row);
				mat.setRow(mat.numRows() - 1, 0, col);
				mat.set(mat.numRows() - 1, mat.numCols() - 1, 0);

				SimpleMatrix rhs = Utils.filled(mat.numRows(), 1, 1);
				for (int i = 0; i < ediisSize - 1; i++) {
					rhs.set(i, Earray[i]);
				}

				double bestE = 0;
				SimpleMatrix bestDIIS = null;
				int n = mat.numRows() - 2;
				for (int i = 0; i <= n; i++) {
					for (int[] tbr : TBRS[i]) {
						try {
							SimpleMatrix newmat =
									removeElementsSquare(mat, tbr);
							SimpleMatrix newrhs =
									removeElementsLinear(rhs, tbr);
							SimpleMatrix tempEdiis =
									addRows(newmat.solve(newrhs), tbr);
							tempEdiis.set(tempEdiis.numRows() - 1, 0);
							boolean nonNegative = !(CommonOps_DDRM.elementMin(
									tempEdiis.getDDRM()) < 0);

							if (nonNegative) {
								double e = finde(tempEdiis);
								if (e < bestE) {
									bestE = e;
									bestDIIS = tempEdiis;
								}
							}

						} catch (SingularMatrixException ignored) {
						}
					}
				}

				SimpleMatrix finalDIIS = bestDIIS;

				SimpleMatrix F = new SimpleMatrix(densityMatrix.numRows(),
						densityMatrix.numCols());

				for (int i = 0; i < finalDIIS.getNumElements() - 1; i++) {
					F = F.plus(Farray[i].scale(finalDIIS.get(i)));
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
				SimpleMatrix mat = new SimpleMatrix(ediisSize, ediisSize);

				for (int i = 0; i < Math.min(Farray.length, numIt + 1); i++) {
					for (int j = i; j < Math.min(Farray.length, numIt + 1);
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

					SimpleMatrix F =
							new SimpleMatrix(densityMatrix.numRows(),
									densityMatrix.numCols());

					SimpleMatrix D =
							new SimpleMatrix(densityMatrix.numRows(),
									densityMatrix.numCols());


					for (int i = 0; i < DIIS.getNumElements() - 1; i++) {
						F = F.plus(Farray[i].scale(DIIS.get(i)));
						D = D.plus(Darray[i].scale(DIIS.get(i)));
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
				} catch (SingularMatrixException e) { // todo fix adrian
					matrices = Utils.symEigen(F);

					E = matrices[1].diag();

					C = matrices[0].transpose();

					double damp = 0.8;
					densityMatrix = calculateDensityMatrix(C).scale(1 - damp)
							.plus(olddensity.scale(damp));
				}
			}

			numIt++;
		}

		findHF();
		findHomo();
		findDipole();

		return this;
	}

	private static SimpleMatrix commutator(SimpleMatrix F, SimpleMatrix D) {
		return F.mult(D).minus(D.mult(F));
	}

	private static SimpleMatrix removeElementsSquare(SimpleMatrix original,
													 int[] indices) {
		SimpleMatrix newarray =
				new SimpleMatrix(original.numRows() - indices.length,
						original.numRows() - indices.length);

		ArrayList<Integer> array = new ArrayList<>();
		for (int i = 0; i < original.numRows(); i++) {
			array.add(i);
		}

		for (int i : indices) {
			array.remove(Integer.valueOf(i));
		}

		int count = 0;

		for (int i : array) {
			int count1 = 0;
			for (int j : array) {
				newarray.set(count, count1, original.get(i, j));
				count1++;
			}

			count++;
		}

		return newarray;
	}

	private static SimpleMatrix removeElementsLinear(SimpleMatrix original,
													 int[] indices) {//get rid
		// of the rows given in indices and return downsized vector

		SimpleMatrix newarray =
				new SimpleMatrix(original.numRows() - indices.length, 1);

		ArrayList<Integer> array = new ArrayList<>();
		for (int i = 0; i < original.numRows(); i++) {
			array.add(i);
		}

		for (int i : indices) {
			array.remove(Integer.valueOf(i));
		}

		int count = 0;

		for (int i : array) {
			newarray.set(count, original.get(i));

			count++;
		}

		return newarray;
	}

	private static SimpleMatrix addRows(SimpleMatrix original,
										int[] indices) { // add zero row at
		// indices

		SimpleMatrix newarray =
				new SimpleMatrix(original.numRows() + indices.length, 1);

		ArrayList<Double> array = new ArrayList<>();

		for (double i : original.getDDRM().data) {
			array.add(i);
		}

		for (int index : indices) {
			array.add(index, 0.0);
		}

		for (int i = 0; i < array.size(); i++) {
			newarray.set(i, array.get(i));
		}

		return newarray;
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

	private double finde(SimpleMatrix ediis) {
		double e = 0;

		for (int a = 0; a < ediis.getNumElements() - 1; a++) {
			e -= Earray[a] * ediis.get(a);
		}

		for (int a = 0; a < ediis.getNumElements() - 1; a++) {
			for (int b = 0; b < ediis.getNumElements() - 1; b++) {
				e += 0.5 * ediis.get(a) * ediis.get(b) * 0.5 *
						B.get(a, b);
			}
		}
		return e;
	}

	private SimpleMatrix calculateDensityMatrix(
			SimpleMatrix c) {//density matrix construction by definition.
		SimpleMatrix densityMatrix = new SimpleMatrix(orbitals.length,
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
				densityMatrix.set(i, j, sum);
			}
		}
		return densityMatrix;
	}

	public SimpleMatrix getE() {
		return E;
	}

	@Override
	public SimpleMatrix alphaDensity() {
		return densityMatrix.scale(0.5);
	}

	@Override
	public SimpleMatrix betaDensity() {
		return densityMatrix.scale(0.5);
	}

	@Override
	public SimpleMatrix densityMatrix() {
		return densityMatrix;
	}
}
