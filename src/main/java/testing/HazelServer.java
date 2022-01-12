package testing;

import com.hazelcast.config.Config;
import com.hazelcast.core.Hazelcast;
import com.hazelcast.core.HazelcastInstance;
import frontend.FrontendConfig;

public class HazelServer {
	public static void main(String[] args) {
		FrontendConfig.init();
		Config config = new Config();
		config.setClusterName("dev");
		config.setProperty("hazelcast.logging.type", "log4j2");

		HazelcastInstance h1 = Hazelcast.newHazelcastInstance(config);
	}
}
