package nddo;

import java.util.Arrays;

public final class NDDOParams {
	private final double[] params; // unfortunately immutable arrays are not a thing in java
	private final double[] aParams;

	public NDDOParams(double alpha, double betas, double betap, double uss,
						 double upp, double zetas, double zetap, double eisol,
						 double gss, double gsp, double hsp, double gpp, double gp2) {
		params = new double[]{alpha, betas, betap, uss, upp, zetas, zetap, eisol, gss, gsp, hsp, gpp, gp2};
		aParams = new double[0];
	}

	/**
	 * Same as the verbose constructor. Clones array passed in so it's essentially pass-by-value.
	 * NOTE: Use this the exact same way you would use the verbose constructor! It's all cloned!
	 * @param params Params array of size 13 by default. More would cause creation of additional parameters array.
	 */
	public NDDOParams(double[] params) {
		this.params = Arrays.copyOfRange(params, 0, 13);
		this.aParams = Arrays.copyOfRange(params, 13, params.length);
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

	// todo make params final
	public void modifyParam(int index, double amnt) {
		params[index] += amnt;
	}

	public double[] toArray() {
		double[] combinedParams = new double[params.length + aParams.length];

		System.arraycopy(params, 0, combinedParams, 0, params.length);
		System.arraycopy(aParams, 0, combinedParams, params.length, aParams.length);

		return combinedParams;
	}

	public double[] getaParams() {
		return aParams.clone(); // todo maybe not a good idea to clone
	}

	@Override
	public NDDOParams clone() {
		return new NDDOParams(toArray());
	} // todo make this a copy constructor instead
}
