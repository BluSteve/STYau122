package frontend;

import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

public class Main {
	public static void main(String[] args) throws IOException {
		Files.createDirectories(Path.of("pastinputs"));
		Files.createDirectories(Path.of("outputs"));

		RunInput input = TxtIO.readInput();
		JsonIO.write(input, "input");

		RunIterator iterator = new RunIterator(input, 100);
		int i = 0;
		for (RunOutput ro : iterator) {
			TxtIO.updateOutput(ro);
			JsonIO.write(ro.getInput(), String.format("pastinputs/%04d-%s", i, ro.getInput().hash));
			JsonIO.write(ro, String.format("outputs/%04d-%s-%s", i, ro.getInput().hash, ro.hash));
			i++;
		}

		System.exit(0);
	}
}
