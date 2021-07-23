package nddoparam;

import nddoparam.mndo.MNDOAtom;
import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.jblas.Solve;
import runcycle.input.RawMolecule;

import java.util.ArrayList;


public class SolutionR extends Solution {

	public double[] integralArray;
	public DoubleMatrix C, F, G, E;
	//H - core matrix, G = 2-electron matrix, F = fock matrix, C = coeffecient
	// matrix (transposed for easier reading), E = eigenvalues
	private DoubleMatrix densityMatrix, B;
	double[] Earray;



	public SolutionR(NDDOAtom[] atoms, int charge) {
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

//		System.out.println( +
//				" Initial diagonalization completed, beginning SCF " +
//				"iterations" +
//				"." +
//				"..");

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
//			System.out.println("numIt = " + numIt);

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


			int min = Math.min(Farray.length + 1, numIt + 2);
			if (commutatorarray[Math.min(Farray.length - 1, numIt)].max() >
					0.01) {
				// if true do EDIIS else DIIS

				DoubleMatrix mat = DoubleMatrix.zeros(min, min);

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

				for (int i = 0; i < Math.min(Farray.length, numIt + 1); i++) {
					rhs.put(i, Earray[i]);
				}

				boolean nonNegative;
				DoubleMatrix tempEdiis;
				DoubleMatrix bestDIIS = null;
				double bestE = 0;

				tempEdiis = Solve.solve(mat, rhs);
				tempEdiis = tempEdiis.put(tempEdiis.rows - 1, 0);
				nonNegative = !(tempEdiis.min() < 0);

				if (nonNegative) {
					double e = finde(tempEdiis);

					bestE = e;
					bestDIIS = tempEdiis.dup();
				}

				for (int i = 0; i < mat.rows - 2; i++) {
					DoubleMatrix newmat =
							removeElementsSquare(mat.dup(), new int[]{i});
					DoubleMatrix newrhs =
							removeElementsLinear(rhs.dup(), new int[]{i});
					tempEdiis = addRows(Solve.solve(newmat, newrhs),
							new int[]{i});
					tempEdiis = tempEdiis.put(tempEdiis.rows - 1, 0);
					nonNegative = !(tempEdiis.min() < 0);

					if (nonNegative) {
						double e = finde(tempEdiis);

						if (e < bestE) {
							bestE = e;
							bestDIIS = tempEdiis;
						}
					}
				}


				for (int i = 0; i < mat.rows - 2; i++) {
					for (int j = i + 1; j < mat.rows - 2; j++) {
						DoubleMatrix newmat =
								removeElementsSquare(mat.dup(),
										new int[]{i, j});
						DoubleMatrix newrhs =
								removeElementsLinear(rhs.dup(),
										new int[]{i, j});
						tempEdiis = addRows(Solve.solve(newmat, newrhs),
								new int[]{i, j});
						tempEdiis = tempEdiis.put(tempEdiis.rows - 1, 0);
						nonNegative = !(tempEdiis.min() < 0);

						if (nonNegative) {
							double e = finde(tempEdiis);

							if (e < bestE) {
								bestE = e;

								bestDIIS = tempEdiis.dup();
							}
						}
					}

				}

				for (int i = 0; i < mat.rows - 2; i++) {

					for (int j = i + 1; j < mat.rows - 2; j++) {

						for (int k = j + 1; k < mat.rows - 2; k++) {
							try {

								DoubleMatrix newmat =
										removeElementsSquare(mat.dup(),
												new int[]{i, j, k});
								DoubleMatrix newrhs =
										removeElementsLinear(rhs.dup(),
												new int[]{i, j, k});
								tempEdiis = addRows(Solve.solve(newmat,
										newrhs),
										new int[]{i, j, k});

								tempEdiis =
										tempEdiis.put(tempEdiis.rows - 1, 0);

								nonNegative = !(tempEdiis.min() < 0);

								if (nonNegative) {
									double e = finde(tempEdiis);

									if (e < bestE) {
										bestE = e;

										bestDIIS = tempEdiis.dup();
									}
								}
							} catch (Exception ex) {
							}
						}
					}

				}

				for (int i = 0; i < mat.rows - 2; i++) {

					for (int j = i + 1; j < mat.rows - 2; j++) {

						for (int k = j + 1; k < mat.rows - 2; k++) {

							for (int l = k + 1; l < mat.rows - 2; l++) {
								try {

									DoubleMatrix newmat =
											removeElementsSquare(mat.dup(),
													new int[]{i, j, k, l});
									DoubleMatrix newrhs =
											removeElementsLinear(rhs.dup(),
													new int[]{i, j, k, l});
									tempEdiis =
											addRows(Solve.solve(newmat,
													newrhs),
													new int[]{i, j, k, l});

									tempEdiis = tempEdiis
											.put(tempEdiis.rows - 1, 0);

									nonNegative = !(tempEdiis.min() < 0);

									if (nonNegative) {
										double e = finde(tempEdiis);

										if (e < bestE) {
											bestE = e;

											bestDIIS = tempEdiis.dup();
										}
									}
								} catch (Exception ex) {
								}
							}
						}
					}

				}


				for (int i = 0; i < mat.rows - 2; i++) {

					for (int j = i + 1; j < mat.rows - 2; j++) {

						for (int k = j + 1; k < mat.rows - 2; k++) {

							for (int l = k + 1; l < mat.rows - 2; l++) {

								for (int m = l + 1; m < mat.rows - 2; m++) {
									try {

										DoubleMatrix newmat =
												removeElementsSquare(mat.dup(),
														new int[]{i, j, k, l,
																m});
										DoubleMatrix newrhs =
												removeElementsLinear(rhs.dup(),
														new int[]{i, j, k, l,
																m});
										tempEdiis = addRows(Solve
														.solve(newmat, newrhs),
												new int[]{i, j, k, l, m});

										tempEdiis = tempEdiis
												.put(tempEdiis.rows - 1, 0);

										nonNegative = !(tempEdiis.min() < 0);

										if (nonNegative) {
											double e =
													finde(tempEdiis);

											if (e < bestE) {
												bestE = e;

												bestDIIS = tempEdiis.dup();
											}
										}
									} catch (Exception ex) {
									}
								}
							}
						}
					}

				}

				for (int i = 0; i < mat.rows - 2; i++) {

					for (int j = i + 1; j < mat.rows - 2; j++) {

						for (int k = j + 1; k < mat.rows - 2; k++) {

							for (int l = k + 1; l < mat.rows - 2; l++) {

								for (int m = l + 1; m < mat.rows - 2; m++) {

									for (int n = m + 1; n < mat.rows - 2; n++) {

										try {

											DoubleMatrix newmat =
													removeElementsSquare(
															mat.dup(),
															new int[]{i, j, k,
																	l, m, n});
											DoubleMatrix newrhs =
													removeElementsLinear(
															rhs.dup(),
															new int[]{i, j, k,
																	l, m, n});
											tempEdiis = addRows(Solve
															.solve(newmat,
																	newrhs),
													new int[]{i, j, k, l, m,
															n});

											tempEdiis = tempEdiis
													.put(tempEdiis.rows - 1,
															0);

											nonNegative =
													!(tempEdiis.min() < 0);

											if (nonNegative) {
												double e = finde(tempEdiis);

												if (e < bestE) {
													bestE = e;

													bestDIIS = tempEdiis.dup();
												}
											}
										} catch (Exception ex) {
										}
									}
								}
							}
						}
					}

				}

				for (int i = 0; i < mat.rows - 2; i++) {

					for (int j = i + 1; j < mat.rows - 2; j++) {

						for (int k = j + 1; k < mat.rows - 2; k++) {

							for (int l = k + 1; l < mat.rows - 2; l++) {

								for (int m = l + 1; m < mat.rows - 2; m++) {

									for (int n = m + 1; n < mat.rows - 2; n++) {

										for (int o = n + 1; o < mat.rows - 2;
											 o++) {
											try {

												DoubleMatrix newmat =
														removeElementsSquare(
																mat.dup(),
																new int[]{i, j,
																		k, l
																		, m,
																		n, o});
												DoubleMatrix newrhs =
														removeElementsLinear(
																rhs.dup(),
																new int[]{i, j,
																		k, l
																		, m,
																		n, o});
												tempEdiis = addRows(Solve
																.solve(newmat,
																		newrhs),
														new int[]{i, j, k, l
																, m,
																n, o});

												tempEdiis = tempEdiis
														.put(tempEdiis.rows - 1,
																0);

												nonNegative =
														!(tempEdiis.min() < 0);

												if (nonNegative) {
													double e =
															finde(tempEdiis);

													if (e < bestE) {
														bestE = e;

														bestDIIS =
																tempEdiis.dup();
													}
												}
											} catch (Exception ex) {
											}
										}
									}
								}
							}
						}
					}
				}
				// todo wtf
				if (bestDIIS == null) bestDIIS = tempEdiis.dup();
				DoubleMatrix finalDIIS = bestDIIS.dup();
//				System.out.println("finalDIIS = " + finalDIIS);
//				double e = 0;
//
//
//				for (int a = 0;
//					 a < DIIS.length - 1;
//					 a++) {
//					e -= Earray[a] *
//							DIIS.get(a);
//				}
//
//				for (int a = 0;
//					 a < DIIS.length - 1;
//					 a++) {
//					for (int b = 0;
//						 b < DIIS.length -
//								 1; b++) {
//						e += 0.5 * DIIS.get(
//								a) *
//								DIIS.get(
//										b) *
//								0.5 *
//								B.get(a,
//										b);
//					}
//				}
//				System.out.println("e = " + e);

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


			}
			else {

				DoubleMatrix mat = DoubleMatrix
						.zeros(min,
								min);

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
//			System.out.println("\n");

		}


//        System.out.println(moleculeName + " SCF completed: " + numIt + "
//        iterations
//        used");

//        System.err.println("Hybrid DIIS took: " + sw.getTime());

		double e = 0;

		for (int j = 0; j < orbitals.length; j++) {
			for (int k = 0; k < orbitals.length; k++) {
				e += 0.5 * densityMatrix.get(j, k) *
						(H.get(j, k) + F.get(j, k));
			}
		}

//        double checksum = 0;
//
//        for (int a = 0; a < atoms.length; a++) {
//            checksum += E(a, index);
//        }
//
//        for (int a = 0; a < atoms.length; a++) {
//            for (int b = a + 1; b < atoms.length; b++) {
//                checksum += E(a, b, index);
//            }
//        }
//
//        if (Math.abs(checksum - e) > 1E-5 || checksum != checksum) {
//            System.err.println ("I knew it!");
//            System.err.println (checksum);
//            System.err.println (e);
//            System.exit(0);
//        }

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

	private double finde(DoubleMatrix tempEdiis) {
		double e = 0;

		for (int a = 0; a < tempEdiis.length - 1; a++) {
			e -= Earray[a] * tempEdiis.get(a);
		}

		for (int a = 0; a < tempEdiis.length - 1; a++) {
			for (int b = 0; b < tempEdiis.length - 1; b++) {
				e += 0.5 * tempEdiis.get(a) * tempEdiis.get(b) * 0.5 *
						B.get(a, b);
			}
		}
		return e;
	}

	private static DoubleMatrix commutator(DoubleMatrix F, DoubleMatrix D) {

		return F.mmul(D).sub(D.mmul(F));
	}

	private static DoubleMatrix removeElementsSquare(DoubleMatrix original,
													 int[] indices) {//remove
		// rows and
		// columns specified in indices and return downsized square matrix
//		System.out.print(Arrays.toString(indices) + ",");
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
		// of the
		// rows given in indices and return downsized vector

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

		for (int i = 0; i < indices.length; i++) {

			array.add(indices[i], 0.0);
		}

		for (int i = 0; i < array.size(); i++) {
			newarray.put(i, array.get(i));
		}

		return newarray;
	}

	@Override
	public SolutionR setRm(RawMolecule rm) {
		this.rm = rm;
		return this;
	}

	public SolutionR clone() {
		NDDOAtom[] newAtoms = new NDDOAtom[atoms.length];
		for (int i = 0; i < atoms.length; i++) {
			newAtoms[i] = new MNDOAtom((MNDOAtom) atoms[i]);
		}
		return new SolutionR(newAtoms, charge);
	}

	private double E(int atomnum, int[][] index) {

		double e = 0;

		for (int i : index[atomnum]) {
			if (i > -1) {
				e += densityMatrix.get(i, i) * orbitals[i].U();
			}
		}

		for (int i : index[atomnum]) {
			for (int j : index[atomnum]) {
				for (int k : index[atomnum]) {
					for (int l : index[atomnum]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {
							e += 0.5 * densityMatrix.get(i, j) *
									(densityMatrix.get(k, l) *
											NDDO6G.OneCenterERI(orbitals[i],
													orbitals[j],
													orbitals[k], orbitals[l])
											- 0.5 * densityMatrix.get(k, l) *
											NDDO6G.OneCenterERI(orbitals[i],
													orbitals[k],
													orbitals[j], orbitals[l]));
						}
					}
				}
			}
		}

		return e;
	}

	private double E(int atomnum1, int atomnum2, int[][] index) {

		double e = 0;

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				if (i != -1 && j != -1) {
					e += densityMatrix.get(i, j) *
							atoms[atomnum2].V(orbitals[i], orbitals[j]);
				}
			}
		}

		for (int k : index[atomnum2]) {
			for (int l : index[atomnum2]) {
				if (k != -1 && l != -1) {
					e += densityMatrix.get(k, l) *
							atoms[atomnum1].V(orbitals[k], orbitals[l]);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int k : index[atomnum2]) {
				if (i != -1 && k != -1) {
					e += 2 * densityMatrix.get(i, k) *
							NDDO6G.beta(orbitals[i], orbitals[k]);
				}
			}
		}

		for (int i : index[atomnum1]) {
			for (int j : index[atomnum1]) {
				for (int k : index[atomnum2]) {
					for (int l : index[atomnum2]) {
						if (i != -1 && j != -1 && k != -1 && l != -1) {
							e += (densityMatrix.get(i, j) *
									densityMatrix.get(k, l) -
									densityMatrix.get(i, k) * 0.5 *
											densityMatrix.get(j, l))
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

	public DoubleMatrix getE() {
		return E;
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
