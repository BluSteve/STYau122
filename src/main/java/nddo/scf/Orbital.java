package nddo.scf;

import nddo.structs.OrbitalProperties;

public abstract class Orbital {
	public OrbitalProperties op;
	public String type;
	public int shell, L, i, j, k;
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
}
