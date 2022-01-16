package testing;

import com.hazelcast.config.Config;
import com.hazelcast.core.Hazelcast;
import frontend.FrontendConfig;
import org.apache.commons.lang3.time.StopWatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.Serializer;
import tools.Utils;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class HazelServer {
	public static double power;

	public static void main(String[] args) {
		FrontendConfig.init();

		power = getPower();

		Config config = new Config();
		config.setClusterName("beJSaHB3AJQVUBa3G7eSptMopuJCRg");
		config.getNetworkConfig().getJoin().getAutoDetectionConfig().setEnabled(false);
		config.getNetworkConfig().getJoin().getMulticastConfig().setEnabled(false);

		Hazelcast.newHazelcastInstance(config);
	}

	private static double getPower() {
		Logger logger = LogManager.getLogger();

		double power;
		try {
			String powerstr = Files.readString(Path.of("power.txt")).strip();
			power = Double.parseDouble(powerstr);
		} catch (IOException e) {
			logger.info("Starting benchmark...");

			RunInput ri = Serializer.gson.fromJson(Utils.getResource("benchmark-molecules.json"), RunInput.class);

			RunIterator iterator = new RunIterator(ri, 1);

			StopWatch sw = StopWatch.createStarted();
			iterator.next();

			power = 1.0 / sw.getTime() * 100000;

			try {
				FileWriter pw = new FileWriter("power.txt");
				pw.write(Double.toString(power));
				pw.close();
			}
			catch (IOException e2) {
				e2.printStackTrace();
			}
		}

		logger.info("The power of this machine is {}", power);
		return power;
	}
}
