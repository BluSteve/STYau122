package testing;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.util.zip.DeflaterOutputStream;
import java.util.zip.InflaterOutputStream;

public class Compressor {
	private static final Logger logger = LogManager.getLogger();

	public static byte[] deflate(String input) {
		logger.info("input length: {}", input.length());

		ByteArrayOutputStream os = new ByteArrayOutputStream();
		try (DeflaterOutputStream dos = new DeflaterOutputStream(os)) {
			dos.write(input.getBytes(StandardCharsets.UTF_8));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		logger.info("deflated length: {}", os.size());

		return os.toByteArray();
	}

	public static String inflate(byte[] deflated) {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		try (OutputStream ios = new InflaterOutputStream(os)) {
			ios.write(deflated);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		return os.toString(StandardCharsets.UTF_8);
	}
}
