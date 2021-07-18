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
	 *  @param g          ParamGradient object that contains the Solution
	 *                   object required.
	 * @param analytical Boolean indicating whether analytical derivatives
	 *                   should be used when available. Should be true by
	 *                   default.
	 */
	public ParamHessianU(ParamGradientU g, boolean analytical) {
		super(g.s, g.datum, g.sExp, analytical);
		this.g = g;
	}

	@Override
	protected ParamGradient constructGPrime(int ZIndex, int paramNum) {
		return new ParamGradientU(
				new SolutionU(Utils.perturbAtomParams(s.atoms,
						s.getUniqueZs()[ZIndex], paramNum), s.charge,
						s.multiplicity),
				datum, (SolutionU) sExp, analytical);
	}
}
