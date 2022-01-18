package frontend;

import com.google.gson.Gson;
import nddo.Config;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.core.config.Configurator;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class FrontendConfig extends Config {
	public static FrontendConfig config;

	public static void init() {
		Gson gson = new Gson();
		FrontendConfig config;
		try {
			config = gson.fromJson(new FileReader("config.json"), FrontendConfig.class);
		} catch (FileNotFoundException e) {
			config = new FrontendConfig();
			LogManager.getLogger().warn("config.json not found, using default config");
		}
		config.setLoggingLevel();

		LogManager.getLogger().info(config);

		nddo.State.config = config;

		FrontendConfig.config = config;
	}

	public int num_runs = 1, starting_run_num = 0;
	public String logging_level = "info";
	public double reconf_power_threshold = 0.1;

	public static void main(String[] args) throws IOException {
		Gson gson = new Gson();
		FileWriter fw = new FileWriter("config.json");
		gson.toJson(new FrontendConfig(), fw);
		fw.close();
	}

	public Config setLoggingLevel() {
		Level level = Level.toLevel(logging_level);
		Configurator.setRootLevel(level);

		return this;
	}

	@Override
	public String toString() {
		Gson gson = new Gson();
		return gson.toJson(this);
	}
}
