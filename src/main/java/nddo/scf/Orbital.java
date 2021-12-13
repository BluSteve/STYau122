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
		this.type = op.getType();
		this.shell = op.getShell();
		this.i = op.geti();
		this.j = op.getj();
		this.k = op.getk();
		this.L = op.getL();
	}
}
