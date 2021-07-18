package nddoparam.param;

import nddoparam.SolutionR;
import scf.Utils;

public class ParamHessianR extends ParamHessian {
	/**
	 * @see ParamHessian
	 */
	public ParamHessianR(SolutionR s, double[] datum, SolutionR sExp,
						 boolean analytical) {
		super(s, datum, sExp, analytical);
		g = new ParamGradientR(s, datum, sExp, analytical);
		g.compute();
	}

	/**
	 * Takes in a pre-computed ParamGradient object.
	 *  @param g          ParamGradient object that contains the Solution
	 *                   object required.
	 * @param analytical Boolean indicating whether analytical derivatives
	 *                   should be used when available. Should be true by
	 *                   default.
	 */
	public ParamHessianR(ParamGradientR g, boolean analytical) {
		super(g.s, g.datum, g.sExp, analytical);
		this.g = g;
	}

	@Override
	protected ParamGradient constructGPrime(int ZIndex, int paramNum) {
		return new ParamGradientR(
				new SolutionR(Utils.perturbAtomParams(s.atoms,
						s.getUniqueZs()[ZIndex], paramNum), s.charge),
				datum, (SolutionR) sExp, analytical);
	}
}
