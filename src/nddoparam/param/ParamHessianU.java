package nddoparam.param;

import nddoparam.SolutionU;
import scf.Utils;

public class ParamHessianU extends ParamHessian {
	/**
	 * @see ParamHessian
	 */
	public ParamHessianU(SolutionU s, double[] datum, SolutionU sExp,
						 boolean analytical) {
		super(s, datum, sExp, analytical);
		g = new ParamGradientU(s, datum, sExp, analytical);
		g.compute();
	}

	/**
	 * Takes in a pre-computed ParamGradient object.
	 * @param g          ParamGradient object that contains the Solution
	 *                   object required.
	 *
	 */
	public ParamHessianU(ParamGradientU g) {
		super(g.s, g.datum, g.sExp, g.analytical);
		this.g = g;
	}

	@Override
	protected ParamGradient constructGPrime(int ZIndex, int paramNum) {
		return new ParamGradientU(
				new SolutionU(Utils.perturbAtomParams(s.atoms,
						s.getRm().mats[ZIndex], paramNum), s.charge,
						s.multiplicity).setRm(s.getRm()),
				datum, (SolutionU) sExp, analytical);
	}
}
