package nddoparam.param;

import nddoparam.SolutionU;
import scf.Utils;

public class ParamHessianU extends ParamHessian {
	public ParamHessianU(SolutionU s,
						 double[] datum, SolutionU sExp,
						 boolean analytical, int[] atomTypes) {
		super(s, datum, sExp, analytical, atomTypes);
		g = new ParamGradientU(s, datum, sExp, analytical, atomTypes);
		g.computeGradients();
	}

	public ParamHessianU(ParamGradientU g, boolean analytical,
						 int[] atomTypes) {
		super(g.s, g.datum, g.sExp, analytical, atomTypes);
		this.g = g;
	}

	@Override
	protected void constructGPrime(int ZIndex, int paramNum) {
		gPrime = new ParamGradientU(
				new SolutionU(Utils.perturbAtomParams(s.atoms,
						s.getUniqueZs()[ZIndex], paramNum), s.charge,
						s.multiplicity),
				datum, (SolutionU) sExp, analytical, atomTypes);
	}
}
