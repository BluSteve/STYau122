package testing;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.zip.DataFormatException;
import java.util.zip.Deflater;
import java.util.zip.Inflater;

public class Compressor {
	private static final Logger logger = LogManager.getLogger();
	public static void main(String[] args) throws DataFormatException, IOException {
		String input=  Files.readString(Path.of("remote-input.json"));


		byte[] inputBytes = input.getBytes(StandardCharsets.UTF_8);
		System.out.println(inputBytes.length);

		int len = 100 + input.length();
		byte[] outputBytes = new byte[len];
		Deflater deflater = new Deflater();
		deflater.setInput(inputBytes);
		deflater.finish();
		int compressedDataLength = deflater.deflate(outputBytes);
		deflater.end();

		Pair<Integer, byte[]> deflate = deflate(input);
		System.out.println("a" + deflate.getRight().length);
		System.out.println("b" + compressedDataLength);

		Inflater inflater = new Inflater();
		inflater.setInput(outputBytes, 0, compressedDataLength);
		byte[] result = new byte[input.length()];
		int resultLength = inflater.inflate(result);
		inflater.end();

		String outputString = new String(result, 0, resultLength, StandardCharsets.UTF_8);

		System.out.println("a" + inflate(deflate.getLeft(), deflate.getRight()).length());
		System.out.println("b" + outputString.length());
		System.out.println();
	}

	public static Pair<Integer, byte[]> deflate(String input) {
		logger.debug("input length: {}", input.length());

		byte[] inputBytes = input.getBytes(StandardCharsets.UTF_8);

		byte[] outputBytes = new byte[1000 + input.length()];
		Deflater deflater = new Deflater();
		deflater.setInput(inputBytes);
		deflater.finish();
		int outputLength = deflater.deflate(outputBytes);
		deflater.end();

		logger.debug("deflated length: {}", outputLength);

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
