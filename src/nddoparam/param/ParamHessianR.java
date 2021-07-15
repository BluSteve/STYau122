package nddoparam.param;

import nddoparam.SolutionR;
import scf.Utils;

public class ParamHessianR extends ParamHessian {
	public ParamHessianR(SolutionR s, String kind, double[] datum,
						 SolutionR sExp, boolean analytical) {
		super(s, kind, datum, sExp, analytical);
		g = new ParamGradientR(s, datum, sExp, analytical);
		g.computeGradients();
	}

	public ParamHessianR(ParamGradientR g, boolean analytical) {
		super(g.s, g.kind, g.datum, g.sExp, analytical);
		this.g = g;
	}

	@Override
	protected void constructGPrime(int ZIndex, int paramNum) {
		gPrime = new ParamGradientR(
				new SolutionR(Utils.perturbAtomParams(s.atoms,
						s.getUniqueZs()[ZIndex], paramNum), s.charge),
				datum,
				(SolutionR) sExp, analytical);
	}
}
