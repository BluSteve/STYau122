package nddo;

import nddo.math.ERI;
import nddo.param.ParamDerivative;
import nddo.scf.GTO;
import nddo.scf.LCGTO;
import nddo.scf.STO6G;
import nddo.structs.OrbitalProperties;

public class NDDO6G extends STO6G implements NDDOOrbital {
	public double beta, U, p0, p1, p2, D1, D2;
	public double gss, gsp, hsp, gpp, gp2, hp2;
	private final NDDOAtom atom;
	private NDDO6G[] orbitalArray;

	public NDDO6G(NDDOAtom atom, OrbitalProperties op, double zeta, double beta, double U) {
		super(op, atom.getCoordinates(), zeta);

		this.atom = atom;
		this.beta = beta;
		this.U = U;
		this.p0 = atom.p0;
		this.p1 = atom.p1;
		this.p2 = atom.p2;
		this.D1 = atom.D1;
		this.D2 = atom.D2;
		this.gss = atom.getParams().getGss();
		this.gsp = atom.getParams().getGsp();
		this.hsp = atom.getParams().getHsp();
		this.gpp = atom.getParams().getGpp();
		this.gp2 = atom.getParams().getGp2();
		this.hp2 = 0.5 * (gpp - gp2);
	}

	private NDDO6G(NDDO6G nddo6G) {
		this(nddo6G.atom, nddo6G.op, nddo6G.zeta, nddo6G.beta, nddo6G.U);
	}

	private NDDO6G(NDDO6G nddo6G, int i, int j, int k) {
		this(nddo6G);

		super.i = i;
		super.j = j;
		super.k = k;
		super.L = i + j + k;
	}

	public NDDO6G(NDDO6G nddo6G, double[] coordinates) {
		this(nddo6G);

		super.coordinates = coordinates;
	}

	public static double beta(NDDO6G a, NDDO6G b) {
		return 0.5 * (a.beta + b.beta) * LCGTO.getS(a, b);
	}

	public static double betaderiv(NDDO6G a, NDDO6G b, int tau) {
		return 0.5 * (a.beta + b.beta) * LCGTO.getSDeriv(a, b, tau);
	}

	public static double betaparamderiv(NDDO6G a, NDDO6G b, int num,
										int type) {
		return 0.5 * (a.beta + b.beta) * ParamDerivative
				.getSderivfinite(a, b, num, type);
	}

	public static double betaderiv2(NDDO6G a, NDDO6G b, int tau1, int tau2) {
		return 0.5 * (a.beta + b.beta) * LCGTO.getSDeriv2(a, b, tau1, tau2);
	}


	public static double getG(NDDO6G a, NDDO6G b, NDDO6G c, NDDO6G d) {
		double[] coeffA = a.decomposition(a.coordinates, c.coordinates);
		double[] coeffB = b.decomposition(a.coordinates, c.coordinates);
		double[] coeffC = c.decomposition(a.coordinates, c.coordinates);
		double[] coeffD = d.decomposition(a.coordinates, c.coordinates);

		double sum2 = 0;

		NDDO6G[] A = a.orbitalArray();
		NDDO6G[] B = b.orbitalArray();
		NDDO6G[] C = c.orbitalArray();
		NDDO6G[] D = d.orbitalArray();

		for (int i = 0; i < coeffA.length; i++)
			for (int j = 0; j < coeffB.length; j++)
				for (int k = 0; k < coeffC.length; k++)
					for (int l = 0; l < coeffD.length; l++)
						if (coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] != 0) {
							sum2 += coeffA[i] * coeffB[j] * coeffC[k] * coeffD[l] *
											ERI.LocalTwoCenterERI(A[i], B[j], C[k], D[l]) * Constants.eV;
						}

		return sum2;
	}

	@Override
	public NDDOAtom getAtom() {
		return this.atom;
	}

	@Override
	public double U() {
		return this.U;
	}

	@Override
	public double beta() {
		return this.beta;
	}

	@Override
	public double[] decomposition(double[] point1, double[] point2) {
		if (this.L == 0) {
			return new double[]{1};
		}
		else if (this.L == 1) {
			double[] zloc = new double[3];
			double val = GTO.R(point1, point2);
			for (int i = 0; i < 3; i++) {
				zloc[i] = (point2[i] - point1[i]) / val;
			}

			double[] yloc = new double[3];
			double scale = Math.sqrt(zloc[0] * zloc[0] + zloc[1] * zloc[1]);

			yloc[0] = -zloc[1] / scale;
			yloc[1] = zloc[0] / scale;
			yloc[2] = 0;

			if (scale == 0) {
				yloc[0] = 0;
				yloc[1] = 1;
			}

			double[] xloc = new double[]{-yloc[1] * zloc[2], yloc[0] * zloc[2],
					zloc[0] * yloc[1] - zloc[1] * yloc[0]};

			if (this.i == 1) {
				return new double[]{xloc[0], yloc[0], zloc[0]};
			}
			else if (this.j == 1) {
				return new double[]{xloc[1], yloc[1], zloc[1]};
			}
			else if (this.k == 1) {
				return new double[]{xloc[2], yloc[2], zloc[2]};
			}
		}
		return null;
	}

	@Override
	public double[] decomposition2(double[] point1, double[] point2) {
		if (this.L == 0) {
			return new double[]{1};
		}

		double x = point2[0] - point1[0];
		double y = point2[1] - point1[1];
		double z = point2[2] - point1[2];

		double Rxz = Math.sqrt(x * x + z * z);

		double R = GTO.R(point1, point2);

		if (this.L == 1) {
			if (this.i == 1) {
				if (Rxz == 0) {
					return new double[]{1, 0, 0};
				}
				return new double[]{x * y / (R * Rxz), -z / Rxz, x / R};
			}
			else if (this.j == 1) {
				return new double[]{-Rxz / R, 0, y / R};
			}
			else {
				if (Rxz == 0) {
					return new double[]{0, 1, 0};
				}
				return new double[]{y * z / (R * Rxz), x / Rxz, z / R};
			}
		}

		return null;
	}

	@Override
	public NDDO6G[] orbitalArray() {
		if (this.orbitalArray == null) {
			if (this.L == 0) {
				orbitalArray = new NDDO6G[]{new NDDO6G(this, 0, 0, 0)};
			}
			else if (this.L == 1) {
				orbitalArray = new NDDO6G[]{
						new NDDO6G(this, 1, 0, 0),
						new NDDO6G(this, 0, 1, 0),
						new NDDO6G(this, 0, 0, 1)};
			}
		}

		return this.orbitalArray;
	}
}
