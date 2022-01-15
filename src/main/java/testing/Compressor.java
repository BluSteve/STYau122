package testing;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.DeflaterOutputStream;
import java.util.zip.InflaterOutputStream;

public class Compressor {
	private static final Logger logger = LogManager.getLogger();

	public static void main(String[] args) throws IOException {
		String json = Files.readString(Path.of("remote-output.json"));

		byte[] d2 = deflate(json);
		System.out.println("deflate: " + d2.length);
		System.out.println(inflate(d2).substring(0, 100));
	}

	public static byte[] deflate(String input) {
		logger.debug("input length: {}", input.length());

		ByteArrayOutputStream os = new ByteArrayOutputStream();
		try (DeflaterOutputStream dos = new DeflaterOutputStream(os)) {
			dos.write(input.getBytes(StandardCharsets.UTF_8));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		logger.debug("deflated length: {}", os.size());

		return os.toByteArray();
	}

	public static String inflate(byte[] deflated) {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		try (OutputStream ios = new InflaterOutputStream(os)) {
			ios.write(deflated);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		logger.debug("inflated length: {}", os.size());

		return os.toString(StandardCharsets.UTF_8);
	}
}
