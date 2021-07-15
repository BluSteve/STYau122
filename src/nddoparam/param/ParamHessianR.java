package nddoparam.param;

import nddoparam.SolutionR;
import scf.Utils;

public class ParamHessianR extends ParamHessian {
	public ParamHessianR(SolutionR s, double[] datum,
						 SolutionR sExp, boolean analytical, int[] atomTypes) {
		super(s, datum, sExp, analytical, atomTypes);
		g = new ParamGradientR(s, datum, sExp, analytical, atomTypes);
		g.computeGradients();
	}

	public ParamHessianR(ParamGradientR g, boolean analytical,
						 int[] atomTypes) {
		super(g.s, g.datum, g.sExp, analytical, atomTypes);
		this.g = g;
	}

	@Override
	protected ParamGradient constructGPrime(int ZIndex, int paramNum) {
		ParamGradient gPrime = new ParamGradientR(
				new SolutionR(Utils.perturbAtomParams(s.atoms,
						s.getUniqueZs()[ZIndex], paramNum), s.charge),
				datum,
				(SolutionR) sExp, analytical, atomTypes);
		return gPrime;
	}
}
