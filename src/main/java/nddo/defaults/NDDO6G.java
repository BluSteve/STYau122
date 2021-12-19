package nddo.defaults;

import nddo.NDDOOrbital;
import nddo.scf.STO6G;
import nddo.structs.OrbitalProperties;

public class NDDO6G extends STO6G implements NDDOOrbital<NDDOAtomBasic, NDDO6G> {
	private final NDDOAtomBasic atom;
	public final double beta, U, p0, p1, p2, D1, D2;
	public final double gss, gsp, hsp, gpp, gp2, hp2;
	private NDDO6G[] orbitalArray;

	public NDDO6G(NDDOAtomBasic atom, OrbitalProperties op, double zeta, double beta, double U) {
		super(op, atom.getCoordinates(), zeta);

		this.atom = atom;
		this.beta = beta;
		this.U = U;
		this.p0 = atom.p0();
		this.p1 = atom.p1();
		this.p2 = atom.p2();
		this.D1 = atom.D1();
		this.D2 = atom.D2();
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

	@Override
	public NDDOAtomBasic getAtom() {
		return this.atom;
	}

	public double U() {
		return this.U;
	}

	@Override
	public NDDO6G[] orbitalArray() { // todo make this better
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
