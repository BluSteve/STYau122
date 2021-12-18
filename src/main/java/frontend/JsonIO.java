package frontend;

import runcycle.structs.Serializer;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


public class JsonIO {
	public static <T> T read(String filename, Class<T> clazz) throws FileNotFoundException {
		return Serializer.gson.fromJson(new FileReader(filename), clazz);
	}

	public static RunInput readInput(String filename) throws FileNotFoundException {
		return Serializer.gson.fromJson(new FileReader(filename), RunInput.class);
	}

	public static RunOutput readOutput(String filename) throws FileNotFoundException {
		return Serializer.gson.fromJson(new FileReader(filename), RunOutput.class);
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
}
