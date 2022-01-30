package remote;

import java.nio.charset.StandardCharsets;

import static runcycle.structs.Serializer.gson;

public class Utils {
	public static byte[] toJsonBytes(Object object) {
		return gson.toJson(object).getBytes(StandardCharsets.UTF_8);
	}

	public static <T> T fromJsonBytes(byte[] bytes, Class<T> clazz) {
		return gson.fromJson(new String(bytes), clazz);
	}
}
