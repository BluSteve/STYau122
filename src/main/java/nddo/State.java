package nddo;

import nddo.geometry.GeometryDerivative;
import nddo.geometry.GeometrySecondDerivative;
import nddo.param.ParamDerivative;

public class State {
	private static ParamDerivative pd = new ParamDerivative();
	private static GeometryDerivative gd = new GeometryDerivative();
	private static GeometrySecondDerivative gd2 = new GeometrySecondDerivative();

	public static ParamDerivative getPd() {
		return pd;
	}

	public static void setPd(ParamDerivative pd) {
		if (pd != null) State.pd = pd;
	}

	public static GeometryDerivative getGd() {
		return gd;
	}

	public static void setGd(GeometryDerivative gd) {
		if (gd != null) State.gd = gd;
	}

	public static GeometrySecondDerivative getGd2() {
		return gd2;
	}

	public static void setGd2(GeometrySecondDerivative gd2) {
		if (gd2 != null) State.gd2 = gd2;
	}
}
