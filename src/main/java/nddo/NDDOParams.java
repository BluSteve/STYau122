package nddo;

public class NDDOParams {
	protected final double[] params; // params are final and read only

	protected NDDOParams(double alpha, double betas, double betap, double uss,
						 double upp, double zetas, double zetap, double eisol,
						 double gss, double gsp, double hsp, double gpp, double gp2) {
		params = new double[]{alpha, betas, betap, uss, upp, zetas, zetap, eisol, gss, gsp, hsp, gpp, gp2};
	}

	/**
	 * Same as the verbose constructor. Clones array passed in so it's essentially pass-by-value.
	 * NOTE: Use this the exact same way you would use the verbose constructor! It's all cloned!
	 * @param params Params array of size 13.
	 */
	protected NDDOParams(double[] params) {
		if (params.length != 13)
			throw new IllegalArgumentException("Invalid number of NDDO params! (" + params.length + ")");

		this.params = params.clone();
	}

	public double getAlpha() {
		return params[0];
	}

	public double getBetas() {
		return params[1];
	}

	public double getBetap() {
		return params[2];
	}

	public double getUss() {
		return params[3];
	}

	public double getUpp() {
		return params[4];
	}

	public double getZetas() {
		return params[5];
	}

	public double getZetap() {
		return params[6];
	}

	public double getEisol() {
		return params[7];
	}

	public double getGss() {
		return params[8];
	}

	public double getGsp() {
		return params[9];
	}

	public double getHsp() {
		return params[10];
	}

	public double getGpp() {
		return params[11];
	}

	public double getGp2() {
		return params[12];
	}

	public void modifyParam(int index, double amnt) {
		params[index] += amnt;
	}

	@Override
	public NDDOParams clone() {
		return new NDDOParams(params);
	}
}
