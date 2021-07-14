package nddoparam.param;

import nddoparam.SolutionU;
import scf.Utils;

public class ParamHessianU extends ParamHessian {
	public ParamHessianU(SolutionU s, String kind,
						 double[] datum, SolutionU sExp,
						 boolean analytical) {
		super(s, kind, datum, sExp, analytical);
		g = new ParamGradientU(s, kind, datum, sExp, analytical);
		g.computeGradients();
	}

	public ParamHessianU(ParamGradientU g, boolean analytical) {
		super(g.s, g.kind, g.datum, g.sExp, analytical);
		this.g = g;
	}

	@Override
	protected void constructGPrime(int ZIndex, int paramNum) {
		gPrime = new ParamGradientU(
				new SolutionU(Utils.perturbAtomParams(s.atoms,
						s.getUniqueZs()[ZIndex], paramNum), s.charge,
						s.multiplicity),
				kind, datum, (SolutionU) sExp, analytical);
	}
}
