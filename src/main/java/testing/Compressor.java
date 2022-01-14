package testing;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

public class Compressor {
	private static final Logger logger = LogManager.getLogger();
	private static final byte[] outputBytes = new byte[500_000_000]; // 500MB input file

	public static Pair<Integer, byte[]> deflate(String input) {
		logger.info("input length: {}", input.length());

		byte[] inputBytes = input.getBytes(StandardCharsets.UTF_8);

		Deflater deflater = new Deflater();
		deflater.setInput(inputBytes);
		deflater.finish();
		int outputLength = deflater.deflate(outputBytes);
		deflater.end();

		logger.info("deflated length: {}", outputLength);

		return Pair.of(input.length(), Arrays.copyOfRange(outputBytes, 0, outputLength));
	}

	public static String inflate(int originalLength, byte[] deflated) {
		try {
			Inflater inflater = new Inflater();
			inflater.setInput(deflated, 0, deflated.length);
			byte[] result = new byte[originalLength];
			int resultLength = inflater.inflate(result);
			inflater.end();

			return new String(result, 0, resultLength, StandardCharsets.UTF_8);
		}
		catch (DataFormatException e) {
			throw new RuntimeException(e);
		}
	}
}
