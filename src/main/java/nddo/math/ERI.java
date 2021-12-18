package nddo.math;

import nddo.NDDO6G;
import nddo.scf.GTO;
import nddo.scf.LCGTO;

import static nddo.math.Multipoles.*;

public class ERI {
	public static double OneCenterERI(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d) {
		if (a.getL() == b.getL() && a.getL() == 0) {/*(ss|??)*/
			if (c.getL() == d.getL() & c.getL() == 0) {/*(gss*/
				return a.gss;
			}
			else if (c.geti() == 1 && d.geti() == 1 || c.getj() == 1 && d.getj() == 1 ||
					c.getk() == 1 && d.getk() == 1) {/*gsp*/
				return a.gsp;
			}
			else return 0;
		}
		else if (a.getL() == 1 && b.getL() == 1) {/*(pp'|??)*/
			if (a.geti() == 1 && b.geti() == 1 || a.getj() == 1 && b.getj() == 1 ||
					a.getk() == 1 && b.getk() == 1) {/*(pp|??)*/
				if (c.getL() == d.getL() && c.getL() == 0) {/*gsp*/
					return a.gsp;
				}
				else if (c.geti() == 1 && d.geti() == 1 || c.getj() == 1 && d.getj() == 1 ||
						c.getk() == 1 && d.getk() == 1)
					if (a.geti() == c.geti() && a.geti() == 1 || a.getj() == c.getj() && a.getj() == 1 ||
							a.getk() == c.getk() && a.getk() == 1) {/*gpp*/
						return a.gpp;
					}
					else {/*gpp'*/
						return a.gp2;
					}
			}
			else if (c.getL() == d.getL() && c.getL() == 1) {/*(pp'|p''p''')*/
				if (Math.abs(LCGTO.getS(a, c) - 1) < 1E-5 && Math.abs(LCGTO.getS(b, d) - 1) < 1E-5) return a.hp2;
				else if (Math.abs(LCGTO.getS(a, d) - 1) < 1E-5 && Math.abs(LCGTO.getS(b, c) - 1) < 1E-5) return a.hp2;
			}
		}
		else if (a.getL() == 0 && c.getL() == 0 &&
				(b.geti() == 1 && d.geti() == 1 || b.getj() == 1 && d.getj() == 1 ||
						b.getk() == 1 && d.getk() == 1)) return a.hsp;
		else if (a.getL() == 0 && d.getL() == 0 &&
				(b.geti() == 1 && c.geti() == 1 || b.getj() == 1 && c.getj() == 1 ||
						b.getk() == 1 && c.getk() == 1)) return a.hsp;
		else if (b.getL() == 0 && c.getL() == 0 &&
				(a.geti() == 1 && d.geti() == 1 || a.getj() == 1 && d.getj() == 1 ||
						a.getk() == 1 && d.getk() == 1)) return a.hsp;
		else if (b.getL() == 0 && d.getL() == 0 &&
				(a.geti() == 1 && c.geti() == 1 || a.getj() == 1 && c.getj() == 1 ||
						a.getk() == 1 && c.getk() == 1)) return a.hsp;
		return 0;
	}

	public static double LocalTwoCenterERI(NDDO6G a, NDDO6G b, NDDO6G c,
										   NDDO6G d) {

		double R = GTO.R(a.getCoords(), c.getCoords());
		//(??|??)
		switch (a.L) {

			case 0://(s?|??)

				switch (b.L) {

					case 0: //(ss|??)

						switch (c.L) {

							case 0: //(ss|s?);

								switch (d.L) {

									case 0://(ss|ss)
										return ssss(a.p0, a.p1, a.p2, a.D1,
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2, R);

									case 1:
										if (d.k == 1) {//(ss|spz)
											return ssspz(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2, R);
										}
										else {//(ss|sppi) = 0
											return 0;
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							case 1: //(ss|p?)
								if (c.k == 1) {//(ss|pz?)

									switch (d.L) {

										case 0://(ss|pzs)
											return ssspz(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2, R);

										case 1:
											if (d.k == 1) {//(ss|pzpz)
												return sspzpz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);
											}
											else {//(ss|pzppi) = 0
												return 0;
											}
										default:
											return 0;
									}
								}
								else {//(ss|ppi?)

									if (d.L == 1 && d.k == 0 && c.i == d.i &&
											c.j == d.j) {//(ss|ppippi)
										return ssppippi(a.p0, a.p1, a.p2, a.D1,
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2, R);
									}
									else {//all others are 0
										return 0;
									}
								}
							default:
								System.err.println("oh no");
								return 0;

						}
					case 1: //(sp|??)

						if (b.k == 1) {//(spz|??)

							switch (c.L) {

								case 0://(spz|s?)

									switch (d.L) {

										case 0://(spz|ss)
											return spzss(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2, R);

										case 1:
											if (d.k == 1) {//(spz|spz)
												return spzspz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.k == 1) {//(spz|pz?)

										switch (d.L) {

											case 0://(spz|pzs)
												return spzspz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);

											case 1:
												if (d.k == 1) {//(spz|pzpz)
													return spzpzpz(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {//(spz|pzppi) = 0
													return 0;
												}
										}
									}
									else {//(spz|ppi?)
										if (d.i == c.i && d.j == c.j &&
												d.k == 0) {
											return spzppippi(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						}
						else {//(sppi|??)

							switch (c.L) {
								case 0://(sppi|s?)
									if (d.i == b.i && d.j == b.j &&
											d.k == 0) {//(sppi|sppi)
										return sppisppi(a.p0, a.p1, a.p2, a.D1,
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2, R);
									}
									else {
										return 0;
									}
								case 1:
									if (c.k == 1) {
										if (d.i == b.i && d.j == b.j &&
												d.k == 0) {//(sppi|pzppi)
											return sppippipz(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									}
									else {
										if (c.i == b.i && c.j == b.j &&
												c.k == 0) {//(sppi|ppi?)
											switch (d.L) {
												case 0:
													return sppisppi(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												case 1:
													if (d.k == 1) {
														return sppippipz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
													else {
														return 0;
													}
												default:
													return 0;
											}
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						}
					default:
						System.err.println("oh no");
						return 0;
				}

			case 1://(p?|??)
				if (a.k == 1) {//(pz?|??)
					switch (b.L) {
						case 0:
							switch (c.L) {

								case 0://(pzs|s?)

									switch (d.L) {

										case 0://(pzs|ss)
											return spzss(a.p0, a.p1, a.p2,
													a.D1,
													a.D2, c.p0, c.p1, c.p2,
													c.D1, c.D2, R);

										case 1:
											if (d.k == 1) {//(pzs|spz)
												return spzspz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.k == 1) {//(pzs|pz?)

										switch (d.L) {

											case 0://(pzs|pzs)
												return spzspz(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);

											case 1:
												if (d.k == 1) {//(pzs|pzpz)
													return spzpzpz(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {//(pzs|pzppi) = 0
													return 0;
												}
										}
									}
									else {//(pzs|ppi?)
										if (d.i == c.i && d.j == c.j &&
												d.k == 0) {
											return spzppippi(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:

							if (b.k == 1) {//(pzpz|??)

								switch (c.L) {

									case 0://(pzpz|s?)

										switch (d.L) {

											case 0://(pzpz|ss)
												return pzpzss(a.p0, a.p1, a.p2,
														a.D1, a.D2, c.p0, c.p1,
														c.p2, c.D1, c.D2, R);

											case 1:
												if (d.k == 1) {//(pzpz|spz)
													return pzpzspz(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {
													return 0;
												}
										}

									case 1:
										if (c.k == 1) {//(pzpz|pz?)

											switch (d.L) {

												case 0://(pzpz|pzs)
													return pzpzspz(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);

												case 1:
													if (d.k == 1) {//(pzpz
														// |pzpz)
														return pzpzpzpz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
													else {//(pzpz|pzppi) = 0
														return 0;
													}
											}
										}
										else {//(pzpz|ppi?)
											if (d.i == c.i && d.j == c.j &&
													d.k == 0) {
												return pzpzppippi(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R);
											}
											else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							}
							else {//(pzppi|??)

								switch (c.L) {
									case 0://(pzppi|s?)
										if (d.i == b.i && d.j == b.j &&
												d.k == 0) {//(pzppi|sppi)
											return ppipzsppi(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									case 1:
										if (c.k == 1) {
											if (d.i == b.i && d.j == b.j &&
													d.k == 0) {//(pzppi|pzppi)
												return ppipzppipz(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R);
											}
											else {
												return 0;
											}
										}
										else {
											if (c.i == b.i && c.j == b.j &&
													c.k == 0) {//(pzppi|ppi?)
												switch (d.L) {
													case 0:
														return ppipzsppi(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													case 1:
														if (d.k == 1) {
															return ppipzppipz(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	R);
														}
														else {
															return 0;
														}
													default:
														return 0;
												}
											}
											else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							}
					}
				}
				else {//(ppi?|??);

					switch (b.L) {
						case 0://(ppis|??)

							switch (c.L) {
								case 0://(ppis|s?)
									if (d.i == a.i && d.j == a.j &&
											d.k == 0) {//(ppis|sppi)
										return sppisppi(a.p0, a.p1, a.p2, a.D1,
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2, R);
									}
									else {
										return 0;
									}
								case 1:
									if (c.k == 1) {
										if (d.i == a.i && d.j == a.j &&
												d.k == 0) {//(ppis|pzppi)
											return sppippipz(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									}
									else {
										if (c.i == a.i && c.j == a.j &&
												c.k == 0) {//(ppis|ppi?)
											switch (d.L) {
												case 0:
													return sppisppi(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												case 1:
													if (d.k == 1) {
														return sppippipz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
													else {
														return 0;
													}
												default:
													return 0;
											}
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:
							if (b.k == 1) {//(ppipz|??)
								switch (c.L) {
									case 0://(ppipz|s?)
										if (d.i == a.i && d.j == a.j &&
												d.k == 0) {//(ppipz|sppi)
											return ppipzsppi(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2, R);
										}
										else {
											return 0;
										}
									case 1:
										if (c.k == 1) {
											if (d.i == a.i && d.j == a.j &&
													d.k == 0) {//(ppipz|pzppi)
												return ppipzppipz(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R);
											}
											else {
												return 0;
											}
										}
										else {
											if (c.i == a.i && c.j == a.j &&
													c.k == 0) {//(ppipz|ppi?)
												switch (d.L) {
													case 0:
														return ppipzsppi(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													case 1:
														if (d.k == 1) {
															return ppipzppipz(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	R);
														}
														else {
															return 0;
														}
													default:
														return 0;
												}
											}
											else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							}
							else {

								switch (c.L) {
									case 0://(ppippi|s?)
										switch (d.L) {
											case 0://(ppippi|ss)
												if (a.i == b.i && a.j == b.j) {
													return ppippiss(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {
													return 0;
												}
											case 1:
												if (d.k == 1 && a.i == b.i &&
														a.j == b.j) {
													return ppippispz(a.p0,
															a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
												else {
													return 0;
												}
										}

									case 1:
										if (c.k == 1) {
											switch (d.L) {
												case 0://(ppippi|pzs)
													if (a.i == b.i &&
															a.j == b.j) {
														return ppippispz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
													else {
														return 0;
													}

												case 1:
													if (d.k == 1 &&
															a.i == b.i && a.j ==
															b.j) {//(ppippi
														// |pzpz)
														return ppippipzpz(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
													else {
														return 0;
													}
											}
										}
										else {
											if (a.i == b.i && a.j ==
													b.j) {//(pxpx|??) or
												// (pypy|??)

												if (c.L == d.L && c.i == d.i &&
														c.j == d.j &&
														c.k == 0) {
													if (a.i == c.i) {
														return ppippippippi(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, R);
													}
													else {
														return pxpxpypy(a.p0,
																a.p1, a.p2,
																a.D1, a.D2,
																c.p0, c.p1,
																c.p2, c.D1,
																c.D2, R);
													}
												}
												else {
													return 0;
												}

											}
											else {//(pxpy|??) or (pypx|??)
												if (c.L == d.L && c.i != d.i &&
														c.j != d.j &&
														c.k == 0) {
													return pxpypxpy(a.p0, a.p1,
															a.p2, a.D1, a.D2,
															c.p0, c.p1, c.p2,
															c.D1, c.D2, R);
												}
											}
										}

								}

							}
					}

				}
		}

		return 0;
	}

	public static double LocalTwoCenterERIpd(NDDO6G a, NDDO6G b,
											 NDDO6G c,
											 NDDO6G d, double D1deriv,
											 double D2deriv,
											 double p1deriv,
											 double p2deriv, int num,
											 int type) {


		double[] A = a.getCoords();
		double[] C = c.getCoords();

		double R = GTO.R(A, C);
		//(??|??)
		switch (a.getL()) {

			case 0://(s?|??)

				switch (b.getL()) {

					case 0: //(ss|??)

						switch (c.getL()) {

							case 0: //(ss|s?);

								switch (d.getL()) {

									case 0://(ss|ss)
										return sssspd(a.p0, a.p1, a.p2,
												a.D1,
												a.D2, c.p0, c.p1, c.p2, c.D1,
												c.D2
												, R, num, D1deriv, D2deriv,
												p1deriv, p2deriv);

									case 1:
										if (d.getk() == 1) {//(ss|spz)
											return ssspzpd(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
										}
										else {//(ss|sppi) = 0
											return 0;
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							case 1: //(ss|p?)
								if (c.getk() == 1) {//(ss|pz?)

									switch (d.getL()) {

										case 0://(ss|pzs)
											return ssspzpd(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);

										case 1:
											if (d.getk() == 1) {//(ss|pzpz)
												return sspzpzpd(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);
											}
											else {//(ss|pzppi) = 0
												return 0;
											}
										default:
											return 0;
									}
								}
								else {//(ss|ppi?)

									if (d.getL() == 1 && d.getk() == 0 &&
											c.geti() == d.geti() &&
											c.getj() == d.getj()) {//(ss
										// |ppippi)
										return ssppippipd(a.p0, a.p1, a.p2,
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2
												, R, num, D1deriv, D2deriv,
												p1deriv, p2deriv);
									}
									else {//all others are 0
										return 0;
									}
								}
							default:
								System.err.println("oh no");
								return 0;

						}
					case 1: //(sp|??)

						if (b.getk() == 1) {//(spz|??)

							switch (c.getL()) {

								case 0://(spz|s?)

									switch (d.getL()) {

										case 0://(spz|ss)
											return spzsspd(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);

										case 1:
											if (d.getk() == 1) {//(spz|spz)
												return spzspzpd(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(spz|pz?)

										switch (d.getL()) {

											case 0://(spz|pzs)
												return spzspzpd(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);

											case 1:
												if (d.getk() == 1) {//(spz
													// |pzpz)
													return spzpzpzpd(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												}
												else {//(spz|pzppi) = 0
													return 0;
												}
										}
									}
									else {//(spz|ppi?)
										if (d.geti() == c.geti() &&
												d.getj() == c.getj() &&
												d.getk() == 0) {
											return spzppippipd(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						}
						else {//(sppi|??)

							switch (c.getL()) {
								case 0://(sppi|s?)
									if (d.geti() == b.geti() &&
											d.getj() == b.getj() &&
											d.getk() == 0) {//(sppi|sppi)
										return sppisppipd(a.p0, a.p1, a.p2,
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2
												, R, num, D1deriv, D2deriv,
												p1deriv, p2deriv);
									}
									else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == b.geti() &&
												d.getj() == b.getj() &&
												d.getk() == 0) {//(sppi|pzppi)
											return sppippipzpd(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
										}
										else {
											return 0;
										}
									}
									else {
										if (c.geti() == b.geti() &&
												c.getj() == b.getj() &&
												c.getk() == 0) {//(sppi|ppi?)
											switch (d.getL()) {
												case 0:
													return sppisppipd(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												case 1:
													if (d.getk() == 1) {
														return sppippipzpd(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													}
													else {
														return 0;
													}
												default:
													return 0;
											}
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						}
					default:
						System.err.println("oh no");
						return 0;
				}

			case 1://(p?|??)
				if (a.getk() == 1) {//(pz?|??)
					switch (b.getL()) {
						case 0:
							switch (c.getL()) {

								case 0://(pzs|s?)

									switch (d.getL()) {

										case 0://(pzs|ss)
											return spzsspd(a.p0, a.p1, a.p2,
													a.D1, a.D2, c.p0, c.p1,
													c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);

										case 1:
											if (d.getk() == 1) {//(pzs|spz)
												return spzspzpd(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);
											}
											else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(pzs|pz?)

										switch (d.getL()) {

											case 0://(pzs|pzs)
												return spzspzpd(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);

											case 1:
												if (d.getk() == 1) {//(pzs
													// |pzpz)
													return spzpzpzpd(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												}
												else {//(pzs|pzppi) = 0
													return 0;
												}
										}
									}
									else {//(pzs|ppi?)
										if (d.geti() == c.geti() &&
												d.getj() == c.getj() &&
												d.getk() == 0) {
											return spzppippipd(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:

							if (b.getk() == 1) {//(pzpz|??)

								switch (c.getL()) {

									case 0://(pzpz|s?)

										switch (d.getL()) {

											case 0://(pzpz|ss)
												return pzpzsspd(a.p0, a.p1,
														a.p2, a.D1, a.D2, c.p0,
														c.p1, c.p2, c.D1, c.D2,
														R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);

											case 1:
												if (d.getk() == 1) {//(pzpz
													// |spz)
													return pzpzspzpd(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												}
												else {
													return 0;
												}
										}

									case 1:
										if (c.getk() == 1) {//(pzpz|pz?)

											switch (d.getL()) {

												case 0://(pzpz|pzs)
													return pzpzspzpd(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
															p2deriv);

												case 1:
													if (d.getk() ==
															1) {//(pzpz|pzpz)
														return pzpzpzpzpd(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													}
													else {//(pzpz|pzppi) = 0
														return 0;
													}
											}
										}
										else {//(pzpz|ppi?)
											if (d.geti() == c.geti() &&
													d.getj() == c.getj() &&
													d.getk() == 0) {
												return pzpzppippipd(a.p0,
														a.p1, a.p2, a.D1, a.D2,
														c.p0, c.p1, c.p2, c.D1,
														c.D2, R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);
											}
											else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							}
							else {//(pzppi|??)

								switch (c.getL()) {
									case 0://(pzppi|s?)
										if (d.geti() == b.geti() &&
												d.getj() == b.getj() &&
												d.getk() == 0) {//(pzppi|sppi)
											return ppipzsppipd(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
										}
										else {
											return 0;
										}
									case 1:
										if (c.getk() == 1) {
											if (d.geti() == b.geti() &&
													d.getj() == b.getj() &&
													d.getk() ==
															0) {//(pzppi|pzppi)
												return ppipzppipzpd(a.p0,
														a.p1, a.p2, a.D1, a.D2,
														c.p0, c.p1, c.p2, c.D1,
														c.D2, R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);
											}
											else {
												return 0;
											}
										}
										else {
											if (c.geti() == b.geti() &&
													c.getj() == b.getj() &&
													c.getk() ==
															0) {//(pzppi|ppi?)
												switch (d.getL()) {
													case 0:
														return ppipzsppipd(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzpd(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	R, num,
																	D1deriv,
																	D2deriv,
																	p1deriv,
																	p2deriv);
														}
														else {
															return 0;
														}
													default:
														return 0;
												}
											}
											else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							}
					}
				}
				else {//(ppi?|??);

					switch (b.getL()) {
						case 0://(ppis|??)

							switch (c.getL()) {
								case 0://(ppis|s?)
									if (d.geti() == a.geti() &&
											d.getj() == a.getj() &&
											d.getk() == 0) {//(ppis|sppi)
										return sppisppipd(a.p0, a.p1, a.p2,
												a.D1, a.D2, c.p0, c.p1, c.p2,
												c.D1, c.D2
												, R, num, D1deriv, D2deriv,
												p1deriv, p2deriv);
									}
									else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == a.geti() &&
												d.getj() == a.getj() &&
												d.getk() == 0) {//(ppis|pzppi)
											return sppippipzpd(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
										}
										else {
											return 0;
										}
									}
									else {
										if (c.geti() == a.geti() &&
												c.getj() == a.getj() &&
												c.getk() == 0) {//(ppis|ppi?)
											switch (d.getL()) {
												case 0:
													return sppisppipd(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												case 1:
													if (d.getk() == 1) {
														return sppippipzpd(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													}
													else {
														return 0;
													}
												default:
													return 0;
											}
										}
										else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:
							if (b.getk() == 1) {//(ppipz|??)
								switch (c.getL()) {
									case 0://(ppipz|s?)
										if (d.geti() == a.geti() &&
												d.getj() == a.getj() &&
												d.getk() == 0) {//(ppipz|sppi)
											return ppipzsppipd(a.p0, a.p1,
													a.p2, a.D1, a.D2, c.p0,
													c.p1, c.p2, c.D1, c.D2
													, R, num, D1deriv, D2deriv,
													p1deriv, p2deriv);
										}
										else {
											return 0;
										}
									case 1:
										if (c.getk() == 1) {
											if (d.geti() == a.geti() &&
													d.getj() == a.getj() &&
													d.getk() ==
															0) {//(ppipz|pzppi)
												return ppipzppipzpd(a.p0,
														a.p1, a.p2, a.D1, a.D2,
														c.p0, c.p1, c.p2, c.D1,
														c.D2, R, num, D1deriv,
														D2deriv, p1deriv,
														p2deriv);
											}
											else {
												return 0;
											}
										}
										else {
											if (c.geti() == a.geti() &&
													c.getj() == a.getj() &&
													c.getk() ==
															0) {//(ppipz|ppi?)
												switch (d.getL()) {
													case 0:
														return ppipzsppipd(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzpd(
																	a.p0, a.p1,
																	a.p2, a.D1,
																	a.D2, c.p0,
																	c.p1, c.p2,
																	c.D1, c.D2,
																	R, num,
																	D1deriv,
																	D2deriv,
																	p1deriv,
																	p2deriv);
														}
														else {
															return 0;
														}
													default:
														return 0;
												}
											}
											else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							}
							else {

								switch (c.getL()) {
									case 0://(ppippi|s?)
										switch (d.getL()) {
											case 0://(ppippi|ss)
												if (a.geti() == b.geti() &&
														a.getj() == b.getj()) {
													return ppippisspd(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												}
												else {
													return 0;
												}
											case 1:
												if (d.getk() == 1 &&
														a.geti() == b.geti() &&
														a.getj() == b.getj()) {
													return ppippispzpd(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												}
												else {
													return 0;
												}
										}

									case 1:
										if (c.getk() == 1) {
											switch (d.getL()) {
												case 0://(ppippi|pzs)
													if (a.geti() == b.geti() &&
															a.getj() ==
																	b.getj()) {
														return ppippispzpd(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													}
													else {
														return 0;
													}

												case 1:
													if (d.getk() == 1 &&
															a.geti() ==
																	b.geti() &&
															a.getj() ==
																	b.getj()) {//(ppippi
														// |pzpz)
														return ppippipzpzpd(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													}
													else {
														return 0;
													}
											}
										}
										else {
											if (a.geti() == b.geti() &&
													a.getj() ==
															b.getj()) {//(pxpx
												// |??) or (pypy|??)

												if (c.getL() == d.getL() &&
														c.geti() == d.geti() &&
														c.getj() == d.getj() &&
														c.getk() == 0) {
													if (a.geti() == c.geti()) {
														return ppippippippipd(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2, R,
																num, D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													}
													else {
														return pxpxpypypd(
																a.p0, a.p1,
																a.p2, a.D1,
																a.D2, c.p0,
																c.p1, c.p2,
																c.D1, c.D2
																, R, num,
																D1deriv,
																D2deriv,
																p1deriv,
																p2deriv);
													}
												}
												else {
													return 0;
												}

											}
											else {//(pxpy|??) or (pypx|??)
												if (c.getL() == d.getL() &&
														c.geti() != d.geti() &&
														c.getj() != d.getj() &&
														c.getk() == 0) {
													return pxpypxpypd(a.p0,
															a.p1, a.p2, a.D1,
															a.D2, c.p0, c.p1,
															c.p2, c.D1, c.D2
															, R, num, D1deriv,
															D2deriv, p1deriv,
															p2deriv);
												}
											}
										}

								}

							}
					}

				}
		}

		return 0;
	}

	public static double LocalTwoCenterERIg2d(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau1, int tau2) {

		double[] A = a.getCoords();
		double[] C = c.getCoords();

		double R = GTO.R(A, C);
		//(??|??)
		switch (a.getL()) {

			case 0://(s?|??)

				switch (b.getL()) {

					case 0: //(ss|??)

						switch (c.getL()) {

							case 0: //(ss|s?);

								switch (d.getL()) {

									case 0://(ss|ss)
										return ssssg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);

									case 1:
										if (d.getk() == 1) {//(ss|spz)
											return ssspzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
										} else {//(ss|sppi) = 0
											return 0;
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							case 1: //(ss|p?)
								if (c.getk() == 1) {//(ss|pz?)

									switch (d.getL()) {

										case 0://(ss|pzs)
											return ssspzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);

										case 1:
											if (d.getk() == 1) {//(ss|pzpz)
												return sspzpzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
											} else {//(ss|pzppi) = 0
												return 0;
											}
										default:
											return 0;
									}
								} else {//(ss|ppi?)

									if (d.getL() == 1 && d.getk() == 0 && c.geti() == d.geti() && c.getj() == d.getj()) {//(ss|ppippi)
										return ssppippig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
									} else {//all others are 0
										return 0;
									}
								}
							default:
								System.err.println("oh no");
								return 0;

						}
					case 1: //(sp|??)

						if (b.getk() == 1) {//(spz|??)

							switch (c.getL()) {

								case 0://(spz|s?)

									switch (d.getL()) {

										case 0://(spz|ss)
											return spzssg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);

										case 1:
											if (d.getk() == 1) {//(spz|spz)
												return spzspzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
											} else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(spz|pz?)

										switch (d.getL()) {

											case 0://(spz|pzs)
												return spzspzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);

											case 1:
												if (d.getk() == 1) {//(spz|pzpz)
													return spzpzpzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
												} else {//(spz|pzppi) = 0
													return 0;
												}
										}
									} else {//(spz|ppi?)
										if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
											return spzppippig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
										} else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						} else {//(sppi|??)

							switch (c.getL()) {
								case 0://(sppi|s?)
									if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(sppi|sppi)
										return sppisppig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
									} else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(sppi|pzppi)
											return sppippipzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
										} else {
											return 0;
										}
									} else {
										if (c.geti() == b.geti() && c.getj() == b.getj() && c.getk() == 0) {//(sppi|ppi?)
											switch (d.getL()) {
												case 0:
													return sppisppig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
												case 1:
													if (d.getk() == 1) {
														return sppippipzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
													} else {
														return 0;
													}
												default:
													return 0;
											}
										} else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						}
					default:
						System.err.println("oh no");
						return 0;
				}

			case 1://(p?|??)
				if (a.getk() == 1) {//(pz?|??)
					switch (b.getL()) {
						case 0:
							switch (c.getL()) {

								case 0://(pzs|s?)

									switch (d.getL()) {

										case 0://(pzs|ss)
											return spzssg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);

										case 1:
											if (d.getk() == 1) {//(pzs|spz)
												return spzspzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
											} else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(pzs|pz?)

										switch (d.getL()) {

											case 0://(pzs|pzs)
												return spzspzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);

											case 1:
												if (d.getk() == 1) {//(pzs|pzpz)
													return spzpzpzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
												} else {//(pzs|pzppi) = 0
													return 0;
												}
										}
									} else {//(pzs|ppi?)
										if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
											return spzppippig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
										} else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:

							if (b.getk() == 1) {//(pzpz|??)

								switch (c.getL()) {

									case 0://(pzpz|s?)

										switch (d.getL()) {

											case 0://(pzpz|ss)
												return pzpzssg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);

											case 1:
												if (d.getk() == 1) {//(pzpz|spz)
													return pzpzspzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
												} else {
													return 0;
												}
										}

									case 1:
										if (c.getk() == 1) {//(pzpz|pz?)

											switch (d.getL()) {

												case 0://(pzpz|pzs)
													return pzpzspzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);

												case 1:
													if (d.getk() == 1) {//(pzpz|pzpz)
														return pzpzpzpzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
													} else {//(pzpz|pzppi) = 0
														return 0;
													}
											}
										} else {//(pzpz|ppi?)
											if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
												return pzpzppippig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
											} else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							} else {//(pzppi|??)

								switch (c.getL()) {
									case 0://(pzppi|s?)
										if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(pzppi|sppi)
											return ppipzsppig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
										} else {
											return 0;
										}
									case 1:
										if (c.getk() == 1) {
											if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(pzppi|pzppi)
												return ppipzppipzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
											} else {
												return 0;
											}
										} else {
											if (c.geti() == b.geti() && c.getj() == b.getj() && c.getk() == 0) {//(pzppi|ppi?)
												switch (d.getL()) {
													case 0:
														return ppipzsppig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
														} else {
															return 0;
														}
													default:
														return 0;
												}
											} else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							}
					}
				} else {//(ppi?|??);

					switch (b.getL()) {
						case 0://(ppis|??)

							switch (c.getL()) {
								case 0://(ppis|s?)
									if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppis|sppi)
										return sppisppig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
									} else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppis|pzppi)
											return sppippipzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
										} else {
											return 0;
										}
									} else {
										if (c.geti() == a.geti() && c.getj() == a.getj() && c.getk() == 0) {//(ppis|ppi?)
											switch (d.getL()) {
												case 0:
													return sppisppig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
												case 1:
													if (d.getk() == 1) {
														return sppippipzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
													} else {
														return 0;
													}
												default:
													return 0;
											}
										} else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:
							if (b.getk() == 1) {//(ppipz|??)
								switch (c.getL()) {
									case 0://(ppipz|s?)
										if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppipz|sppi)
											return ppipzsppig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
										} else {
											return 0;
										}
									case 1:
										if (c.getk() == 1) {
											if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppipz|pzppi)
												return ppipzppipzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
											} else {
												return 0;
											}
										} else {
											if (c.geti() == a.geti() && c.getj() == a.getj() && c.getk() == 0) {//(ppipz|ppi?)
												switch (d.getL()) {
													case 0:
														return ppipzsppig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
														} else {
															return 0;
														}
													default:
														return 0;
												}
											} else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							} else {

								switch (c.getL()) {
									case 0://(ppippi|s?)
										switch (d.getL()) {
											case 0://(ppippi|ss)
												if (a.geti() == b.geti() && a.getj() == b.getj()) {
													return ppippissg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
												} else {
													return 0;
												}
											case 1:
												if (d.getk() == 1 && a.geti() == b.geti() && a.getj() == b.getj()) {
													return ppippispzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
												} else {
													return 0;
												}
										}

									case 1:
										if (c.getk() == 1) {
											switch (d.getL()) {
												case 0://(ppippi|pzs)
													if (a.geti() == b.geti() && a.getj() == b.getj()) {
														return ppippispzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
													} else {
														return 0;
													}

												case 1:
													if (d.getk() == 1 && a.geti() == b.geti() && a.getj() == b.getj()) {//(ppippi|pzpz)
														return ppippipzpzg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
													} else {
														return 0;
													}
											}
										} else {
											if (a.geti() == b.geti() && a.getj() == b.getj()) {//(pxpx|??) or (pypy|??)

												if (c.getL() == d.getL() && c.geti() == d.geti() && c.getj() == d.getj() && c.getk() == 0) {
													if (a.geti() == c.geti()) {
														return ppippippippig2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
													} else {
														return pxpxpypyg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
													}
												} else {
													return 0;
												}

											} else {//(pxpy|??) or (pypx|??)
												if (c.getL() == d.getL() && c.geti() != d.geti() && c.getj() != d.getj() && c.getk() == 0) {
													return pxpypxpyg2d(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau1, tau2);
												}
											}
										}

								}

							}
					}

				}
		}

		return 0;
	}

	public static double LocalTwoCenterERIgd(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d, int tau) {

		double[] A = a.getCoords();
		double[] C = c.getCoords();

		double R = GTO.R(A, C);

		//(??|??)
		switch (a.getL()) {

			case 0://(s?|??)

				switch (b.getL()) {

					case 0: //(ss|??)

						switch (c.getL()) {

							case 0: //(ss|s?);

								switch (d.getL()) {

									case 0://(ss|ss)
										return ssssgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);

									case 1:
										if (d.getk() == 1) {//(ss|spz)
											return ssspzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
										} else {//(ss|sppi) = 0
											return 0;
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							case 1: //(ss|p?)
								if (c.getk() == 1) {//(ss|pz?)

									switch (d.getL()) {

										case 0://(ss|pzs)
											return ssspzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);

										case 1:
											if (d.getk() == 1) {//(ss|pzpz)
												return sspzpzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
											} else {//(ss|pzppi) = 0
												return 0;
											}
										default:
											return 0;
									}
								} else {//(ss|ppi?)

									if (d.getL() == 1 && d.getk() == 0 && c.geti() == d.geti() && c.getj() == d.getj()) {//(ss|ppippi)
										return ssppippigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
									} else {//all others are 0
										return 0;
									}
								}
							default:
								System.err.println("oh no");
								return 0;

						}
					case 1: //(sp|??)

						if (b.getk() == 1) {//(spz|??)

							switch (c.getL()) {

								case 0://(spz|s?)

									switch (d.getL()) {

										case 0://(spz|ss)
											return spzssgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);

										case 1:
											if (d.getk() == 1) {//(spz|spz)
												return spzspzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
											} else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(spz|pz?)

										switch (d.getL()) {

											case 0://(spz|pzs)
												return spzspzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);

											case 1:
												if (d.getk() == 1) {//(spz|pzpz)
													return spzpzpzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
												} else {//(spz|pzppi) = 0
													return 0;
												}
										}
									} else {//(spz|ppi?)
										if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
											return spzppippigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
										} else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						} else {//(sppi|??)

							switch (c.getL()) {
								case 0://(sppi|s?)
									if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(sppi|sppi)
										return sppisppigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
									} else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(sppi|pzppi)
											return sppippipzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
										} else {
											return 0;
										}
									} else {
										if (c.geti() == b.geti() && c.getj() == b.getj() && c.getk() == 0) {//(sppi|ppi?)
											switch (d.getL()) {
												case 0:
													return sppisppigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
												case 1:
													if (d.getk() == 1) {
														return sppippipzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
													} else {
														return 0;
													}
												default:
													return 0;
											}
										} else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						}
					default:
						System.err.println("oh no");
						return 0;
				}

			case 1://(p?|??)
				if (a.getk() == 1) {//(pz?|??)
					switch (b.getL()) {
						case 0:
							switch (c.getL()) {

								case 0://(pzs|s?)

									switch (d.getL()) {

										case 0://(pzs|ss)
											return spzssgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);

										case 1:
											if (d.getk() == 1) {//(pzs|spz)
												return spzspzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
											} else {
												return 0;
											}
									}

								case 1:
									if (c.getk() == 1) {//(pzs|pz?)

										switch (d.getL()) {

											case 0://(pzs|pzs)
												return spzspzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);

											case 1:
												if (d.getk() == 1) {//(pzs|pzpz)
													return spzpzpzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
												} else {//(pzs|pzppi) = 0
													return 0;
												}
										}
									} else {//(pzs|ppi?)
										if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
											return spzppippigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
										} else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:

							if (b.getk() == 1) {//(pzpz|??)

								switch (c.getL()) {

									case 0://(pzpz|s?)

										switch (d.getL()) {

											case 0://(pzpz|ss)
												return pzpzssgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);

											case 1:
												if (d.getk() == 1) {//(pzpz|spz)
													return pzpzspzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
												} else {
													return 0;
												}
										}

									case 1:
										if (c.getk() == 1) {//(pzpz|pz?)

											switch (d.getL()) {

												case 0://(pzpz|pzs)
													return pzpzspzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);

												case 1:
													if (d.getk() == 1) {//(pzpz|pzpz)
														return pzpzpzpzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
													} else {//(pzpz|pzppi) = 0
														return 0;
													}
											}
										} else {//(pzpz|ppi?)
											if (d.geti() == c.geti() && d.getj() == c.getj() && d.getk() == 0) {
												return pzpzppippigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
											} else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							} else {//(pzppi|??)

								switch (c.getL()) {
									case 0://(pzppi|s?)
										if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(pzppi|sppi)
											return ppipzsppigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
										} else {
											return 0;
										}
									case 1:
										if (c.getk() == 1) {
											if (d.geti() == b.geti() && d.getj() == b.getj() && d.getk() == 0) {//(pzppi|pzppi)
												return ppipzppipzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
											} else {
												return 0;
											}
										} else {
											if (c.geti() == b.geti() && c.getj() == b.getj() && c.getk() == 0) {//(pzppi|ppi?)
												switch (d.getL()) {
													case 0:
														return ppipzsppigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
														} else {
															return 0;
														}
													default:
														return 0;
												}
											} else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}
							}
					}
				} else {//(ppi?|??);

					switch (b.getL()) {
						case 0://(ppis|??)

							switch (c.getL()) {
								case 0://(ppis|s?)
									if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppis|sppi)
										return sppisppigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
									} else {
										return 0;
									}
								case 1:
									if (c.getk() == 1) {
										if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppis|pzppi)
											return sppippipzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
										} else {
											return 0;
										}
									} else {
										if (c.geti() == a.geti() && c.getj() == a.getj() && c.getk() == 0) {//(ppis|ppi?)
											switch (d.getL()) {
												case 0:
													return sppisppigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
												case 1:
													if (d.getk() == 1) {
														return sppippipzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
													} else {
														return 0;
													}
												default:
													return 0;
											}
										} else {
											return 0;
										}
									}
								default:
									System.err.println("oh no");
									return 0;
							}
						case 1:
							if (b.getk() == 1) {//(ppipz|??)
								switch (c.getL()) {
									case 0://(ppipz|s?)
										if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppipz|sppi)
											return ppipzsppigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
										} else {
											return 0;
										}
									case 1:
										if (c.getk() == 1) {
											if (d.geti() == a.geti() && d.getj() == a.getj() && d.getk() == 0) {//(ppipz|pzppi)
												return ppipzppipzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
											} else {
												return 0;
											}
										} else {
											if (c.geti() == a.geti() && c.getj() == a.getj() && c.getk() == 0) {//(ppipz|ppi?)
												switch (d.getL()) {
													case 0:
														return ppipzsppigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
													case 1:
														if (d.getk() == 1) {
															return ppipzppipzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
														} else {
															return 0;
														}
													default:
														return 0;
												}
											} else {
												return 0;
											}
										}
									default:
										System.err.println("oh no");
										return 0;
								}

							} else {

								switch (c.getL()) {
									case 0://(ppippi|s?)
										switch (d.getL()) {
											case 0://(ppippi|ss)
												if (a.geti() == b.geti() && a.getj() == b.getj()) {
													return ppippissgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
												} else {
													return 0;
												}
											case 1:
												if (d.getk() == 1 && a.geti() == b.geti() && a.getj() == b.getj()) {
													return ppippispzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
												} else {
													return 0;
												}
										}

									case 1:
										if (c.getk() == 1) {
											switch (d.getL()) {
												case 0://(ppippi|pzs)
													if (a.geti() == b.geti() && a.getj() == b.getj()) {
														return ppippispzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
													} else {
														return 0;
													}

												case 1:
													if (d.getk() == 1 && a.geti() == b.geti() && a.getj() == b.getj()) {//(ppippi|pzpz)
														return ppippipzpzgd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
													} else {
														return 0;
													}
											}
										} else {
											if (a.geti() == b.geti() && a.getj() == b.getj()) {//(pxpx|??) or (pypy|??)

												if (c.getL() == d.getL() && c.geti() == d.geti() && c.getj() == d.getj() && c.getk() == 0) {
													if (a.geti() == c.geti()) {
														return ppippippippigd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
													} else {
														return pxpxpypygd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
													}
												} else {
													return 0;
												}

											} else {//(pxpy|??) or (pypx|??)
												if (c.getL() == d.getL() && c.geti() != d.geti() && c.getj() != d.getj() && c.getk() == 0) {
													return pxpypxpygd(a.p0, a.p1, a.p2, a.D1, a.D2, c.p0, c.p1, c.p2, c.D1, c.D2, R, A, C, tau);
												}
											}
										}

								}

							}
					}

				}
		}

		return 0;
	}
}
