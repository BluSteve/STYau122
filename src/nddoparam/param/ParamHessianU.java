package nddoparam.param;

import nddoparam.SolutionU;
import scf.Utils;

public class ParamHessianU extends ParamHessian {
	/**
	 * @see ParamHessian
	 */
	public ParamHessianU(SolutionU s, double[] datum, SolutionU sExp,
						 boolean analytical, int[] atomTypes) {
		super(s, datum, sExp, analytical);
		g = new ParamGradientU(s, datum, sExp, analytical);
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
	public ParamHessianU(ParamGradientU g, boolean analytical,
						 int[] atomTypes) {
		super(g.s, g.datum, g.sExp, analytical);
		this.g = g;
	}

	@Override
	protected ParamGradient constructGPrime(int ZIndex, int paramNum) {
		ParamGradient gPrime = new ParamGradientU(
				new SolutionU(Utils.perturbAtomParams(s.atoms,
						s.getUniqueZs()[ZIndex], paramNum), s.charge,
						s.multiplicity),
				datum, (SolutionU) sExp, analytical);
		return gPrime;
	}
}
