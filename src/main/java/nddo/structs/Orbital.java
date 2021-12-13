package nddo.structs;

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
		this.L = op.getL();
		this.i = op.getConfig()[0];
		this.j = op.getConfig()[1];
		this.k = op.getConfig()[2];
	}
}
