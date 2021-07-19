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
	 *
	 */
	public ParamHessianR(ParamGradientR g) {
		super(g.s, g.datum, g.sExp, g.analytical);
		this.g = g;
	}

	@Override
	protected ParamGradient constructGPrime(int ZIndex, int paramNum) {
		return new ParamGradientR(
				new SolutionR(Utils.perturbAtomParams(s.atoms,
						s.getRm().mats[ZIndex], paramNum), s.charge)
						.setRm(s.getRm()),
				datum, (SolutionR) sExp, analytical);
	}
}
