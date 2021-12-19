package nddo;

public final class NDDOParams { // params finality is up to the user
	public final double[] params;

	/**
	 * Same as the verbose constructor. Clones array passed in so it's essentially pass-by-value.
	 * NOTE: Use this the exact same way you would use the verbose constructor! It's all cloned!
	 * @param params Params array.
	 */
	public NDDOParams(double ...params) {
		if (params == null) throw new NullPointerException("Params cannot be null!");
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

	public double get(int index) {
		return params[index];
	}

	/**
	 * Gets index corresponding to additional params.
	 * @param index i
	 * @return params[13+i]
	 */
	public double aget(int index) {
		return params[index + 13];
	}

	public NDDOParams copy() {
		return new NDDOParams(params);
	}
}
