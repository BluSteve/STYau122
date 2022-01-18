package frontend;

import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;
import runcycle.structs.Serializer;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class JsonIO {
	private static final ExecutorService es = Executors.newSingleThreadExecutor();

	public static <T> T read(String filename, Class<T> clazz) throws FileNotFoundException {
		return Serializer.gson.fromJson(new FileReader(filename+ ".json"), clazz);
	}

	public static RunInput readInput(String filename) throws FileNotFoundException {
		return Serializer.gson.fromJson(new FileReader(filename + ".json"), RunInput.class);
	}

	public static RunOutput readOutput(String filename) throws FileNotFoundException {
		return Serializer.gson.fromJson(new FileReader(filename+ ".json"), RunOutput.class);
	}

	public static void write(Object o, String... filenames) throws IOException {
		if (filenames.length == 0) {
			FileWriter fw = new FileWriter(Serializer.getHash(o) + ".json");
			Serializer.gson.toJson(o, fw);
			fw.close();
		}
		else for (String filename : filenames) {
			FileWriter fw = new FileWriter(filename + ".json");
			Serializer.gson.toJson(o, fw);
			fw.close();
		}
	}

	public static void writeAsync(Object o, String... filenames) {
		es.submit(() -> {
			try {
				write(o, filenames);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		});
	}
}
