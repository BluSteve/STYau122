//package examples.am1;
//
//import nddo.NDDOParams;
//
//import java.util.Arrays;
//
//public class AM1Params  {
//	protected final double[] params2;
//
//	public AM1Params(double alpha, double betas, double betap, double uss, double upp, double zetas, double zetap,
//					 double eisol, double gss, double gsp, double hsp, double gpp, double gp2,
//					 double K1, double K2, double K3, double K4, double L1, double L2, double L3, double L4,
//					 double M1, double M2, double M3, double M4) {
//		super(alpha, betas, betap, uss, upp, zetas, zetap, eisol, gss, gsp, hsp, gpp, gp2);
//		params2 = new double[]{K1, K2, K3, K4, L1, L2, L3, L4, M1, M2, M3, M4};
//	}
//
//	/**
//	 * Initializes AM1Params from an array, copying all elements.
//	 *
//	 * @param params 25-element array.
//	 */
//	public AM1Params(double[] params) {
//		super(Arrays.copyOfRange(params, 0, 13));
//
//		if (params.length != 25)
//			throw new IllegalArgumentException("Invalid number of AM1 params! (" + params.length + ")");
//
//		params2 = Arrays.copyOfRange(params, 13, params.length);
//	}
//
//	public double getK1() {
//		return params2[0];
//	}
//
//	public double getK2() {
//		return params2[1];
//	}
//
//	public double getK3() {
//		return params2[2];
//	}
//
//	public double getK4() {
//		return params2[3];
//	}
//
//	public double getL1() {
//		return params2[4];
//	}
//
//	public double getL2() {
//		return params2[5];
//	}
//
//	public double getL3() {
//		return params2[6];
//	}
//
//	public double getL4() {
//		return params2[7];
//	}
//
//	public double getM1() {
//		return params2[8];
//	}
//
//	public double getM2() {
//		return params2[9];
//	}
//
//	public double getM3() {
//		return params2[10];
//	}
//
//	public double getM4() {
//		return params2[11];
//	}
//
//	@Override
//	public AM1Params clone() {
//		double[] combinedParams = new double[params.length + params2.length];
//
//		System.arraycopy(params, 0, combinedParams, 0, params.length);
//		System.arraycopy(params2, 0, combinedParams, params.length, params2.length);
//
//		return new AM1Params(combinedParams);
//	}
//}