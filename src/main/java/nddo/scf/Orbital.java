package nddo.scf;

import nddo.structs.OrbitalProperties;

public abstract class Orbital {
	public final OrbitalProperties op;
	public final String type;
	public final int shell;
	public int L, i, j, k;
	public double[] coordinates;

	protected Orbital(OrbitalProperties op, double[] coordinates) {
		this.coordinates = coordinates;
		this.op = op;
		this.type = op.type;
		this.shell = op.shell;
		this.i = op.i;
		this.j = op.j;
		this.k = op.k;
		this.L = op.L;
	}

	protected Orbital(Orbital orbital) { // not a copy constructor!
		this.coordinates = orbital.coordinates;
		this.op = orbital.op;
		this.type = op.type;
		this.shell = op.shell;
		this.i = op.i;
		this.j = op.j;
		this.k = op.k;
		this.L = op.L;
	}

	public double[] derivativeDecomposition(double[] point1, double[] point2, int tau) {
		if (L == 0) return new double[]{0};
		double R = GTO.R(point1, point2);
		double Rxy = Math.sqrt(
				(point2[1] - point1[1]) * (point2[1] - point1[1]) + (point2[0] - point1[0]) * (point2[0] - point1[0]));
		switch (tau) {
			case 0:
				if (i == 1) {
					double x1 = (point2[2] - point1[2]) / (R * Rxy) -
							(point2[0] - point1[0]) * (point2[0] - point1[0]) * (point2[2] - point1[2]) /
									(Rxy * Rxy * Rxy * R) -
							(point2[0] - point1[0]) * (point2[0] - point1[0]) * (point2[2] - point1[2]) /
									(R * R * R * Rxy);
					double x2 = -(point2[0] - point1[0]) * (point2[1] - point1[1]) / (Rxy * Rxy * Rxy);
					double x3 = (point2[0] - point1[0]) * (point2[0] - point1[0]) / (R * R * R) - 1 / R;
					return new double[]{x1, x2, x3};
				}
				else if (j == 1) {
					double x1 = -(point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2]) /
							(Rxy * Rxy * Rxy * R) -
							(point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2]) /
									(R * R * R * Rxy);
					double x2 = (point2[0] - point1[0]) * (point2[0] - point1[0]) / (Rxy * Rxy * Rxy) - 1 / Rxy;
					double x3 = (point2[0] - point1[0]) * (point2[1] - point1[1]) / (R * R * R);
					return new double[]{x1, x2, x3};
				}
				else if (k == 1) {
					double x1 = (point2[0] - point1[0]) * Rxy / (R * R * R) - (point2[0] - point1[0]) / (R * Rxy);
					double x2 = 0;
					double x3 = (point2[0] - point1[0]) * (point2[2] - point1[2]) / (R * R * R);
					return new double[]{x1, x2, x3};
				}
			case 1:
				if (i == 1) {
					double x1 = -(point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2]) /
							(Rxy * Rxy * Rxy * R) -
							(point2[0] - point1[0]) * (point2[1] - point1[1]) * (point2[2] - point1[2]) /
									(R * R * R * Rxy);
					double x2 = -(point2[1] - point1[1]) * (point2[1] - point1[1]) / (Rxy * Rxy * Rxy) + 1 / Rxy;
					double x3 = (point2[0] - point1[0]) * (point2[1] - point1[1]) / (R * R * R);
					return new double[]{x1, x2, x3};
				}
				else if (j == 1) {
					double x1 = (point2[2] - point1[2]) / (R * Rxy) -
							(point2[2] - point1[2]) * (point2[1] - point1[1]) * (point2[1] - point1[1]) /
									(Rxy * Rxy * Rxy * R) -
							(point2[2] - point1[2]) * (point2[1] - point1[1]) * (point2[1] - point1[1]) /
									(R * R * R * Rxy);
					double x2 = (point2[0] - point1[0]) * (point2[1] - point1[1]) / (Rxy * Rxy * Rxy);
					double x3 = (point2[1] - point1[1]) * (point2[1] - point1[1]) / (R * R * R) - 1 / R;
					return new double[]{x1, x2, x3};
				}
				else if (k == 1) {
					double x1 = (point2[1] - point1[1]) * Rxy / (R * R * R) - (point2[1] - point1[1]) / (R * Rxy);
					double x2 = 0;
					double x3 = (point2[2] - point1[2]) * (point2[1] - point1[1]) / (R * R * R);
					return new double[]{x1, x2, x3};
				}
			case 2:
				if (i == 1) {
					double x1 = (point2[0] - point1[0]) / (R * Rxy) -
							(point2[0] - point1[0]) * (point2[2] - point1[2]) * (point2[2] - point1[2]) /
									(R * R * R * Rxy);
					double x2 = 0;
					double x3 = (point2[2] - point1[2]) * (point2[0] - point1[0]) / (R * R * R);
					return new double[]{x1, x2, x3};
				}
				else if (j == 1) {
					double x1 = (point2[1] - point1[1]) / (R * Rxy) -
							(point2[1] - point1[1]) * (point2[2] - point1[2]) * (point2[2] - point1[2]) /
									(R * R * R * Rxy);
					double x2 = 0;
					double x3 = (point2[2] - point1[2]) * (point2[1] - point1[1]) / (R * R * R);
					return new double[]{x1, x2, x3};
				}
				else if (k == 1) {
					double x1 = (point2[2] - point1[2]) * Rxy / (R * R * R);
					double x2 = 0;
					double x3 = (point2[2] - point1[2]) * (point2[2] - point1[2]) / (R * R * R) - 1 / R;
					return new double[]{x1, x2, x3};
				}
		}
		return null;
	}

	public double[] derivativeDecompositionvar(double[] point1, double[] point2, int tau) {
		if (L == 0) return new double[]{0};
		double x = point2[0] - point1[0];
		double y = point2[1] - point1[1];
		double z = point2[2] - point1[2];
		double Rxz = Math.sqrt(x * x + z * z);
		double R = GTO.R(point1, point2);
		if (L == 1) switch (tau) {
			case 0:
				if (i == 1) {
					double val1 = x * x * y / (R * R * R * Rxz) + x * x * y / (R * Rxz * Rxz * Rxz) - y / (R * Rxz);
					double val2 = -x * z / (Rxz * Rxz * Rxz);
					double val3 = x * x / (R * R * R) - 1 / R;
					return new double[]{val1, val2, val3};
				}
				else if (j == 1) {
					double val1 = x / (R * Rxz) - x * Rxz / (R * R * R);
					double val3 = x * y / (R * R * R);
					return new double[]{val1, 0, val3};
				}
				else if (k == 1) {
					double val1 = x * y * z / (R * R * R * Rxz) + x * y * z / (R * Rxz * Rxz * Rxz);
					double val2 = x * x / (Rxz * Rxz * Rxz) - 1 / Rxz;
					double val3 = x * z / (R * R * R);
					return new double[]{val1, val2, val3};
				}
			case 1:
				if (i == 1) {
					double val1 = x * y * y / (R * R * R * Rxz) - x / (R * Rxz);
					double val3 = x * y / (R * R * R);
					return new double[]{val1, 0, val3};
				}
				else if (j == 1) {
					double val1 = -y * Rxz / (R * R * R);
					double val3 = y * y / (R * R * R) - 1 / R;
					return new double[]{val1, 0, val3};
				}
				else if (k == 1) {
					double val1 = y * y * z / (R * R * R * Rxz) - z / (R * Rxz);
					double val3 = y * z / (R * R * R);
					return new double[]{val1, 0, val3};
				}
			case 2:
				if (i == 1) {
					double val1 = x * y * z / (R * R * R * Rxz) + x * y * z / (R * Rxz * Rxz * Rxz);
					double val2 = 1 / Rxz - z * z / (Rxz * Rxz * Rxz);
					double val3 = x * z / (R * R * R);
					return new double[]{val1, val2, val3};
				}
				else if (j == 1) {
					double val1 = z / (R * Rxz) - z * Rxz / (R * R * R);
					double val3 = z * y / (R * R * R);
					return new double[]{val1, 0, val3};
				}
				else if (k == 1) {
					double val1 = y * z * z / (R * R * R * Rxz) + y * z * z / (R * Rxz * Rxz * Rxz) - y / (R * Rxz);
					double val2 = x * z / (Rxz * Rxz * Rxz);
					double val3 = z * z / (R * R * R) - 1 / R;
					return new double[]{val1, val2, val3};
				}
		}
		return null;
	}

	public double[] secondDerivativeDecomposition(double[] point1, double[] point2, int tau1, int tau2) {
		if (L == 0) return new double[]{0};
		int A = Math.min(tau1, tau2);
		int B = Math.max(tau1, tau2);
		double x = point2[0] - point1[0];
		double y = point2[1] - point1[1];
		double z = point2[2] - point1[2];
		double R = GTO.R(point1, point2);
		double Rxy = Math.sqrt(x * x + y * y);
		double[] returnval = new double[3];
		if (i == 1 && L == 1) switch (A) {
			case 0:
				switch (B) {
					case 0:/*partial wrt x and x*/
						returnval[0] = 3 * x * z / (R * R * R * Rxy) + 3 * x * z / (R * Rxy * Rxy * Rxy) -
								2 * x * x * x * z / (R * R * R * Rxy * Rxy * Rxy) -
								3 * x * x * x * z / (R * Rxy * Rxy * Rxy * Rxy * Rxy) -
								3 * x * x * x * z / (R * R * R * R * R * Rxy);
						returnval[1] = y / (Rxy * Rxy * Rxy) - 3 * x * x * y / (Rxy * Rxy * Rxy * Rxy * Rxy);
						returnval[2] = 3 * x * x * x / (R * R * R * R * R) - 3 * x / (R * R * R);
						return returnval;
					case 1:/*partial wrt x and y*/
						returnval[0] = z * y / (R * R * R * Rxy) + z * y / (R * Rxy * Rxy * Rxy) -
								2 * x * x * y * z / (R * R * R * Rxy * Rxy * Rxy) -
								3 * x * x * y * z / (R * R * R * R * R * Rxy) -
								3 * x * x * y * z / (Rxy * Rxy * Rxy * Rxy * Rxy * R);
						returnval[1] = x / (Rxy * Rxy * Rxy) - 3 * x * y * y / (Rxy * Rxy * Rxy * Rxy * Rxy);
						returnval[2] = 3 * x * x * y / (R * R * R * R * R) - y / (R * R * R);
						return returnval;
					case 2:/*partial wrt x and z*/
						returnval[0] = z * z / (R * R * R * Rxy) - 1 / (R * Rxy) + x * x / (R * R * R * Rxy) +
								x * x / (Rxy * Rxy * Rxy * R) - 3 * x * x * z * z / (R * R * R * R * R * Rxy) -
								x * x * z * z / (R * R * R * Rxy * Rxy * Rxy);
						returnval[1] = 0;
						returnval[2] = 3 * x * x * z / (R * R * R * R * R) - z / (R * R * R);
						return returnval;
				}
				break;
			case 1:
				switch (B) {
					case 1:/*partial wrt y and y*/
						returnval[0] = x * z / (R * R * R * Rxy) + x * z / (R * Rxy * Rxy * Rxy) -
								2 * x * y * y * z / (R * R * R * Rxy * Rxy * Rxy) -
								3 * x * y * y * z / (R * R * R * R * R * Rxy) -
								3 * x * y * y * z / (Rxy * Rxy * Rxy * Rxy * Rxy * R);
						returnval[1] = 3 * y / (Rxy * Rxy * Rxy) - 3 * y * y * y / (Rxy * Rxy * Rxy * Rxy * Rxy);
						returnval[2] = 3 * x * y * y / (R * R * R * R * R) - x / (R * R * R);
						return returnval;
					case 2:/*partial wrt y and z*/
						returnval[0] = x * y / (R * R * R * Rxy) + x * y / (R * Rxy * Rxy * Rxy) -
								3 * x * y * z * z / (R * R * R * R * R * Rxy) -
								x * y * z * z / (R * R * R * Rxy * Rxy * Rxy);
						returnval[1] = 0;
						returnval[2] = 3 * x * y * z / (R * R * R * R * R);
						return returnval;
				}
			case 2:
				returnval[0] = 3 * x * z / (R * R * R * Rxy) - 3 * x * z * z * z / (R * R * R * R * R * Rxy);
				returnval[1] = 0;
				returnval[2] = 3 * x * z * z / (R * R * R * R * R) - x / (R * R * R);
				return returnval;
			default:
		}
		else if (j == 1 && L == 1) switch (A) {
			case 0:
				switch (B) {
					case 0: /*partial wrt x and x*/
						returnval[0] = y * z / (R * R * R * Rxy) + y * z / (R * Rxy * Rxy * Rxy) -
								2 * x * x * y * z / (R * R * R * Rxy * Rxy * Rxy) -
								3 * x * x * y * z / (R * R * R * R * R * Rxy) -
								3 * x * x * y * z / (Rxy * Rxy * Rxy * Rxy * Rxy * R);
						returnval[1] = -3 * x / (Rxy * Rxy * Rxy) + 3 * x * x * x / (Rxy * Rxy * Rxy * Rxy * Rxy);
						returnval[2] = 3 * x * x * y / (R * R * R * R * R) - y / (R * R * R);
						return returnval;
					case 1:/*partial wrt x and y*/
						returnval[0] = z * x / (R * R * R * Rxy) + z * x / (R * Rxy * Rxy * Rxy) -
								2 * x * y * y * z / (R * R * R * Rxy * Rxy * Rxy) -
								3 * x * y * y * z / (R * R * R * R * R * Rxy) -
								3 * x * y * y * z / (Rxy * Rxy * Rxy * Rxy * Rxy * R);
						returnval[1] = -y / (Rxy * Rxy * Rxy) + 3 * x * x * y / (Rxy * Rxy * Rxy * Rxy * Rxy);
						returnval[2] = 3 * x * y * y / (R * R * R * R * R) - x / (R * R * R);
						return returnval;
					case 2:/*partial wrt x and z*/
						returnval[0] = x * y / (R * R * R * Rxy) + x * y / (R * Rxy * Rxy * Rxy) -
								3 * x * y * z * z / (R * R * R * R * R * Rxy) -
								x * y * z * z / (R * R * R * Rxy * Rxy * Rxy);
						returnval[1] = 0;
						returnval[2] = 3 * x * y * z / (R * R * R * R * R);
						return returnval;
				}
				break;
			case 1:
				switch (B) {
					case 1: /*partial wrt y and y*/
						returnval[0] = 3 * y * z / (R * R * R * Rxy) + 3 * y * z / (R * Rxy * Rxy * Rxy) -
								2 * y * y * y * z / (R * R * R * Rxy * Rxy * Rxy) -
								3 * y * y * y * z / (R * Rxy * Rxy * Rxy * Rxy * Rxy) -
								3 * y * y * y * z / (R * R * R * R * R * Rxy);
						returnval[1] = -x / (Rxy * Rxy * Rxy) + 3 * x * y * y / (Rxy * Rxy * Rxy * Rxy * Rxy);
						returnval[2] = 3 * y * y * y / (R * R * R * R * R) - 3 * y / (R * R * R);
						return returnval;
					case 2:/*partial wrt y and z*/
						returnval[0] = z * z / (R * R * R * Rxy) - 1 / (R * Rxy) + y * y / (R * R * R * Rxy) +
								y * y / (Rxy * Rxy * Rxy * R) - 3 * y * y * z * z / (R * R * R * R * R * Rxy) -
								y * y * z * z / (R * R * R * Rxy * Rxy * Rxy);
						returnval[1] = 0;
						returnval[2] = 3 * y * y * z / (R * R * R * R * R) - z / (R * R * R);
						return returnval;
				}
			case 2:
				returnval[0] = 3 * y * z / (R * R * R * Rxy) - 3 * y * z * z * z / (R * R * R * R * R * Rxy);
				returnval[1] = 0;
				returnval[2] = 3 * y * z * z / (R * R * R * R * R) - y / (R * R * R);
				return returnval;
			default:
		}
		else if (k == 1 && L == 1) switch (A) {
			case 0:
				switch (B) {
					case 0: /*partial wrt x and x*/
						returnval[0] = 3 * Rxy * x * x / (R * R * R * R * R) - Rxy / (R * R * R) + 1 / (Rxy * R) -
								2 * x * x / (Rxy * R * R * R) - x * x / (R * Rxy * Rxy * Rxy);
						returnval[1] = 0;
						returnval[2] = 3 * x * x * z / (R * R * R * R * R) - z / (R * R * R);
						return returnval;
					case 1:/*partial wrt x and y*/
						returnval[0] = 3 * Rxy * x * y / (R * R * R * R * R) - 2 * x * y / (Rxy * R * R * R) -
								x * y / (R * Rxy * Rxy * Rxy);
						returnval[1] = 0;
						returnval[2] = 3 * x * y * z / (R * R * R * R * R);
						return returnval;
					case 2:/*partial wrt x and z*/
						returnval[0] = 3 * Rxy * x * z / (R * R * R * R * R) - x * z / (Rxy * R * R * R);
						returnval[1] = 0;
						returnval[2] = 3 * x * z * z / (R * R * R * R * R) - x / (R * R * R);
						return returnval;
				}
				break;
			case 1:
				switch (B) {
					case 1: /*partial wrt y and y*/
						returnval[0] = 3 * Rxy * y * y / (R * R * R * R * R) - Rxy / (R * R * R) + 1 / (Rxy * R) -
								2 * y * y / (Rxy * R * R * R) - y * y / (R * Rxy * Rxy * Rxy);
						returnval[1] = 0;
						returnval[2] = 3 * y * y * z / (R * R * R * R * R) - z / (R * R * R);
						return returnval;
					case 2:/*partial wrt y and z*/
						returnval[0] = 3 * Rxy * y * z / (R * R * R * R * R) - y * z / (Rxy * R * R * R);
						returnval[1] = 0;
						returnval[2] = 3 * y * z * z / (R * R * R * R * R) - y / (R * R * R);
						return returnval;
				}
			case 2:
				returnval[0] = 3 * z * z * Rxy / (R * R * R * R * R) - Rxy / (R * R * R);
				returnval[1] = 0;
				returnval[2] = 3 * z * z * z / (R * R * R * R * R) - 3 * z / (R * R * R);
				return returnval;
			default:
		}
		System.err.println("oh no!");
		return new double[]{0, 0, 0};
	}

	public double[] secondDerivativeDecompositionvar(double[] point1, double[] point2, int tau1, int tau2) {
		if (L == 0) return new double[]{0};
		int A = Math.min(tau1, tau2);
		int B = Math.max(tau1, tau2);
		double x = point2[0] - point1[0];
		double y = point2[1] - point1[1];
		double z = point2[2] - point1[2];
		double R = GTO.R(point1, point2);
		double Rxz = Math.sqrt(x * x + z * z);
		double[] returnval = new double[3];
		if (i == 1 && L == 1) switch (A) {
			case 0:
				switch (B) {
					case 0:/*partial wrt x and x*/
						returnval[0] = 3 * x * x * x * y / (R * Rxz * Rxz * Rxz * Rxz * Rxz) +
								3 * x * x * x * y / (R * R * R * R * R * Rxz) +
								2 * x * x * x * y / (R * R * R * Rxz * Rxz * Rxz) - 3 * x * y / (R * R * R * Rxz) -
								3 * x * y / (R * Rxz * Rxz * Rxz);
						returnval[1] = z / (Rxz * Rxz * Rxz) - 3 * x * x * z / (Rxz * Rxz * Rxz * Rxz * Rxz);
						returnval[2] = 3 * x * x * x / (R * R * R * R * R) - 3 * x / (R * R * R);
						return returnval;
					case 1:/*partial wrt x and y*/
						returnval[0] = 1 / (R * Rxz) + 3 * x * x * y * y / (R * R * R * R * R * Rxz) +
								x * x * y * y / (R * R * R * Rxz * Rxz * Rxz) - x * x / (R * Rxz * Rxz * Rxz) -
								x * x / (R * R * R * Rxz) - y * y / (R * R * R * Rxz);
						returnval[1] = 0;
						returnval[2] = 3 * x * x * y / (R * R * R * R * R) - y / (R * R * R);
						return returnval;
					case 2:/*partial wrt x and z*/
						returnval[0] = 3 * x * x * y * z / (R * R * R * R * R * Rxz) +
								2 * x * x * y * z / (R * R * R * Rxz * Rxz * Rxz) +
								3 * x * x * y * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz) - y * z / (R * R * R * Rxz) -
								y * z / (R * Rxz * Rxz * Rxz);
						returnval[1] = x / (Rxz * Rxz * Rxz) - 3 * x * z * z / (Rxz * Rxz * Rxz * Rxz * Rxz);
						returnval[2] = 3 * x * x * z / (R * R * R * R * R) - z / (R * R * R);
						return returnval;
				}
				break;
			case 1:
				switch (B) {
					case 1:/*partial wrt y and y*/
						returnval[0] = 3 * x * y * y * y / (R * R * R * R * R * Rxz) - 3 * x * y / (R * R * R * Rxz);
						returnval[1] = 0;
						returnval[2] = 3 * x * y * y / (R * R * R * R * R) - x / (R * R * R);
						return returnval;
					case 2:/*partial wrt y and z*/
						returnval[0] = 3 * x * y * y * z / (R * R * R * R * R * Rxz) +
								x * y * y * z / (R * R * R * Rxz * Rxz * Rxz) - x * z / (R * Rxz * Rxz * Rxz) -
								x * z / (R * R * R * Rxz);
						returnval[1] = 0;
						returnval[2] = 3 * x * y * z / (R * R * R * R * R);
						return returnval;
				}
			case 2:
				returnval[0] = 3 * x * y * z * z / (R * R * R * R * R * Rxz) +
						3 * x * y * z * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz) +
						2 * x * y * z * z / (R * R * R * Rxz * Rxz * Rxz) - x * y / (R * Rxz * Rxz * Rxz) -
						x * y / (R * R * R * Rxz);
				returnval[1] = 3 * z / (Rxz * Rxz * Rxz) - 3 * z * z * z / (Rxz * Rxz * Rxz * Rxz * Rxz);
				returnval[2] = 3 * x * z * z / (R * R * R * R * R) - x / (R * R * R);
				return returnval;
			default:
		}
		else if (j == 1 && L == 1) switch (A) {
			case 0:
				switch (B) {
					case 0: /*partial wrt x and x*/
						returnval[0] =
								2 * x * x / (R * R * R * Rxz) + x * x / (R * Rxz * Rxz * Rxz) + Rxz / (R * R * R) -
										1 / (R * Rxz) - 3 * x * x * Rxz / (R * R * R * R * R);
						returnval[1] = 0;
						returnval[2] = 3 * x * x * y / (R * R * R * R * R) - y / (R * R * R);
						return returnval;
					case 1:/*partial wrt x and y*/
						returnval[0] = x * y / (R * R * R * Rxz) - 3 * x * y * Rxz / (R * R * R * R * R);
						returnval[1] = 0;
						returnval[2] = 3 * x * y * y / (R * R * R * R * R) - x / (R * R * R);
						return returnval;
					case 2:/*partial wrt x and z*/
						returnval[0] = 2 * x * z / (R * R * R * Rxz) + x * z / (R * Rxz * Rxz * Rxz) -
								3 * x * z * Rxz / (R * R * R * R * R);
						returnval[1] = 0;
						returnval[2] = 3 * x * y * z / (R * R * R * R * R);
						return returnval;
				}
				break;
			case 1:
				switch (B) {
					case 1: /*partial wrt y and y*/
						returnval[0] = Rxz / (R * R * R) - 3 * y * y * Rxz / (R * R * R * R * R);
						returnval[1] = 0;
						returnval[2] = 3 * y * y * y / (R * R * R * R * R) - 3 * y / (R * R * R);
						return returnval;
					case 2:/*partial wrt y and z*/
						returnval[0] = y * z / (R * R * R * Rxz) - 3 * y * z * Rxz / (R * R * R * R * R);
						returnval[1] = 0;
						returnval[2] = 3 * y * y * z / (R * R * R * R * R) - z / (R * R * R);
						return returnval;
				}
			case 2:
				returnval[0] = 2 * z * z / (R * R * R * Rxz) + z * z / (R * Rxz * Rxz * Rxz) + Rxz / (R * R * R) -
						1 / (R * Rxz) - 3 * z * z * Rxz / (R * R * R * R * R);
				returnval[1] = 0;
				returnval[2] = 3 * y * z * z / (R * R * R * R * R) - y / (R * R * R);
				return returnval;
			default:
		}
		else if (k == 1 && L == 1) switch (A) {
			case 0:
				switch (B) {
					case 0: /*partial wrt x and x*/
						returnval[0] = 2 * x * x * y * z / (R * R * R * Rxz * Rxz * Rxz) +
								3 * x * x * y * z / (R * R * R * R * R * Rxz) +
								3 * x * x * y * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz) - y * z / (R * R * R * Rxz) -
								y * z / (R * Rxz * Rxz * Rxz);
						returnval[1] = 3 * x * x * x / (Rxz * Rxz * Rxz * Rxz * Rxz) - 3 * x / (Rxz * Rxz * Rxz);
						returnval[2] = 3 * x * x * z / (R * R * R * R * R) - z / (R * R * R);
						return returnval;
					case 1:/*partial wrt x and y*/
						returnval[0] = 3 * x * y * y * z / (R * R * R * R * R * Rxz) +
								x * y * y * z / (R * R * R * Rxz * Rxz * Rxz) - x * z / (R * Rxz * Rxz * Rxz) -
								x * z / (R * R * R * Rxz);
						returnval[1] = 0;
						returnval[2] = 3 * x * y * z / (R * R * R * R * R);
						return returnval;
					case 2:/*partial wrt x and z*/
						returnval[0] = 2 * x * y * z * z / (R * R * R * Rxz * Rxz * Rxz) +
								3 * x * y * z * z / (R * R * R * R * R * Rxz) +
								3 * x * y * z * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz) - x * y / (R * Rxz * Rxz * Rxz) -
								x * y / (R * R * R * Rxz);
						returnval[1] = 3 * x * x * z / (Rxz * Rxz * Rxz * Rxz * Rxz) - z / (Rxz * Rxz * Rxz);
						returnval[2] = 3 * x * z * z / (R * R * R * R * R) - x / (R * R * R);
						return returnval;
				}
				break;
			case 1:
				switch (B) {
					case 1: /*partial wrt y and y*/
						returnval[0] = 3 * y * y * y * z / (R * R * R * R * R * Rxz) - 3 * y * z / (R * R * R * Rxz);
						returnval[1] = 0;
						returnval[2] = 3 * y * y * z / (R * R * R * R * R) - z / (R * R * R);
						return returnval;
					case 2:/*partial wrt y and z*/
						returnval[0] = 1 / (R * Rxz) + 3 * y * y * z * z / (R * R * R * R * R * Rxz) +
								y * y * z * z / (R * R * R * Rxz * Rxz * Rxz) - y * y / (R * R * R * Rxz) -
								z * z / (R * Rxz * Rxz * Rxz) - z * z / (R * R * R * Rxz);
						returnval[1] = 0;
						returnval[2] = 3 * y * z * z / (R * R * R * R * R) - y / (R * R * R);
						return returnval;
				}
			case 2:
				returnval[0] = 3 * y * z * z * z / (R * R * R * R * R * Rxz) +
						2 * y * z * z * z / (R * R * R * Rxz * Rxz * Rxz) +
						3 * y * z * z * z / (R * Rxz * Rxz * Rxz * Rxz * Rxz) - 3 * y * z / (R * Rxz * Rxz * Rxz) -
						3 * y * z / (R * R * R * Rxz);
				returnval[1] = 3 * x * z * z / (Rxz * Rxz * Rxz * Rxz * Rxz) - x / (Rxz * Rxz * Rxz);
				returnval[2] = 3 * z * z * z / (R * R * R * R * R) - 3 * z / (R * R * R);
				return returnval;
			default:
		}
		System.err.println("oh no!");
		return new double[]{0, 0, 0};
	}

	public double[] decomposition(double[] point1, double[] point2) {
		if (this.L == 0) return new double[]{1};
		else if (this.L == 1) {
			double[] zloc = new double[3];
			double val = GTO.R(point1, point2);
			for (int i = 0; i < 3; i++) zloc[i] = (point2[i] - point1[i]) / val;
			double[] yloc = new double[3];
			double scale = Math.sqrt(zloc[0] * zloc[0] + zloc[1] * zloc[1]);
			yloc[0] = -zloc[1] / scale;
			yloc[1] = zloc[0] / scale;
			yloc[2] = 0;
			if (scale == 0) {
				yloc[0] = 0;
				yloc[1] = 1;
			}
			double[] xloc = new double[]{-yloc[1] * zloc[2], yloc[0] * zloc[2], zloc[0] * yloc[1] - zloc[1] * yloc[0]};
			if (this.i == 1) return new double[]{xloc[0], yloc[0], zloc[0]};
			else if (this.j == 1) return new double[]{xloc[1], yloc[1], zloc[1]};
			else if (this.k == 1) return new double[]{xloc[2], yloc[2], zloc[2]};
		}
		return null;
	}

	public double[] decompositionvar(double[] point1, double[] point2) {
		if (this.L == 0) return new double[]{1};
		double x = point2[0] - point1[0];
		double y = point2[1] - point1[1];
		double z = point2[2] - point1[2];
		double Rxz = Math.sqrt(x * x + z * z);
		double R = GTO.R(point1, point2);
		if (this.L == 1) if (this.i == 1) {
			if (Rxz == 0) return new double[]{1, 0, 0};
			return new double[]{x * y / (R * Rxz), -z / Rxz, x / R};
		}
		else if (this.j == 1) return new double[]{-Rxz / R, 0, y / R};
		else {
			if (Rxz == 0) return new double[]{0, 1, 0};
			return new double[]{y * z / (R * Rxz), x / Rxz, z / R};
		}
		return null;
	}
}
