package nddo;

import com.google.gson.Gson;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

import java.io.FileWriter;
import java.io.IOException;

public class Config {
	public int rhf_numIt_max = 200;
	public int uhf_numIt_max = 400;

	public double poplethiel_tolerable = 1e-5, poplethiel_ideal = 1e-10;
	public int poplethiel_batch_size = 1;
	public int rhf_thiel_limit = 10000, uhf_thiel_limit = 10000;
	public int geom_opt_hessian_interval = 7;

	public String logging_level = "info";

	public static void main(String[] args) throws IOException {
		Gson gson = new Gson();
		FileWriter fw = new FileWriter("config.json");
		gson.toJson(new Config(), fw);
		fw.close();
	}

	public void setLoggingLevel() {
		Level level = Level.toLevel(logging_level);
		Configurator.setRootLevel(level);
	}

	@Override
	public String toString() {
		Gson gson = new Gson();
		return gson.toJson(this);
	}
}
