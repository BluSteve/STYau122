package nddo;

import com.google.gson.Gson;
import nddo.defaults.NDDO6GMethods;
import nddo.geometry.GeometryDerivative;
import nddo.geometry.GeometrySecondDerivative;
import nddo.param.ParamDerivative;
import org.apache.logging.log4j.LogManager;

import java.io.FileNotFoundException;
import java.io.FileReader;

public class State {
	public static NDDOOrbitalMethods nom = new NDDO6GMethods();
	public static ParamDerivative pd = new ParamDerivative();
	public static GeometryDerivative gd = new GeometryDerivative();
	public static GeometrySecondDerivative gd2 = new GeometrySecondDerivative();
	public static final Config config = getConfig();

	private static Config getConfig() {
		Gson gson = new Gson();
		Config config;
		try {
			config = gson.fromJson(new FileReader("config.json"), Config.class);
		} catch (FileNotFoundException e) {
			config = new Config();
			LogManager.getLogger().warn("config.json not found, using default config");
		}
		config.setLoggingLevel();

		LogManager.getLogger().info(config);
		return config;
	}
}
