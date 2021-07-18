package nddoparam.param;

import nddoparam.SolutionR;
import scf.Utils;

public class ParamHessianR extends ParamHessian {
	/**
	 * @see ParamHessian
	 */
	public ParamHessianR(SolutionR s, double[] datum, SolutionR sExp,
						 boolean analytical, int[] atomTypes) {
		super(s, datum, sExp, analytical);
		g = new ParamGradientR(s, datum, sExp, analytical);
		g.compute();
	}

	/**
	 * Takes in a pre-computed ParamGradient object.
	 *
	 * @param g          ParamGradient object that contains the Solution
	 *                   object required.
	 * @param analytical Whether to use analytical derivatives if possible.
	 * @param analytical Boolean indicating whether analytical derivatives
	 *                   should be used when available. Should be true by
	 *                   default.
	 * @param atomTypes  The atom types contained in the training set. This is
	 *                   required on construction due to the Hessian matrix
	 *                   having to be constructed with padding for speed
	 *                   enhancements.
	 */
	public ParamHessianR(ParamGradientR g, boolean analytical,
						 int[] atomTypes) {
		super(g.s, g.datum, g.sExp, analytical);
		this.g = g;
	}

	@Override
	protected ParamGradient constructGPrime(int ZIndex, int paramNum) {
		ParamGradient gPrime = new ParamGradientR(
				new SolutionR(Utils.perturbAtomParams(s.atoms,
						s.getUniqueZs()[ZIndex], paramNum), s.charge),
				datum, (SolutionR) sExp, analytical);
		return gPrime;
	}
}
