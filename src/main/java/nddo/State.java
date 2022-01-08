package nddo;

import nddo.defaults.NDDO6GMethods;
import nddo.geometry.GeometryDerivative;
import nddo.geometry.GeometrySecondDerivative;
import nddo.param.ParamDerivative;

public class State {
	public static NDDOOrbitalMethods nom = new NDDO6GMethods();
	public static ParamDerivative pd = new ParamDerivative();
	public static GeometryDerivative gd = new GeometryDerivative();
	public static GeometrySecondDerivative gd2 = new GeometrySecondDerivative();
	public static Config config = new Config();
}
