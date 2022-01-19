package analysis;

import runcycle.structs.RunOutput;
import runcycle.structs.Serializer;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.stream.IntStream;

public class Analysis {
	public static void main(String[] args) throws IOException {
		File dir = new File("outputsgradient/outputs");
		File[] files = dir.listFiles();
		double[] errors = new double[300];
		IntStream.range(0, 300).parallel().forEach( i-> {
			System.out.println(i);
			RunOutput ro = null;
			try {
				ro = Serializer.gson.fromJson(Files.readString(files[i].toPath()), RunOutput.class);
			} catch (IOException e) {
				e.printStackTrace();
			}
			errors[i] = ro.ttError;
		});

		for (double error : errors) {
			System.out.println(error + ",");
		}
	}
}
