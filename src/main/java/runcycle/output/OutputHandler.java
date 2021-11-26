package runcycle.output;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonSyntaxException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import runcycle.MoleculeResult;
import runcycle.input.RawInput;
import tools.Utils;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class OutputHandler {
	private static final Logger logger = LogManager.getLogger();

	public static MoleculeOutput toMoleculeOutput(MoleculeResult result,
												  boolean isRunHessian) {
		MoleculeOutput mo = new MoleculeOutput();
		mo.rawMolecule = result.getRm();
		mo.time = result.getTime();

		mo.hf = result.getHF();
		mo.dipole = result.getDipole();
		mo.ie = -result.getIE();
		mo.geomGradient = result.getGeomGradient();
		mo.totalError = result.getTotalError();

		ParamGradientOutput pgo = new ParamGradientOutput();
		pgo.hf = result.getHFDerivs();
		pgo.dipole = result.getDipoleDerivs();
		pgo.ie = result.getIEDerivs();
		pgo.geom = result.getGeomDerivs();
		pgo.total = result.getTotalGradients();

		mo.gradient = pgo;

		if (isRunHessian) mo.hessian = result.getHessian();

		return mo;
	}

	public static MoleculeOutput[] importMoleculeOutputs(String inputPath) {
		MoleculeOutput[] mos = null;
		try {
			mos = new Gson().fromJson(new FileReader(inputPath + ".json"),
					MoleculeOutput[].class);
		} catch (JsonSyntaxException | FileNotFoundException ignored) {
		}
		return mos;
	}

	// assumes updated input
	public static void output(RawInput ri, MoleculeOutput[] mos,
							  String output, long ttTime) throws IOException {
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		double ttError = 0;
		StringBuilder hashsb = new StringBuilder();
		for (MoleculeOutput mo : mos) {
			ttError += mo.totalError;

			// remove time to make hashes consistent
			long t = mo.time;
			mo.time = 0;
			hashsb.append(gson.toJson(mo));
			mo.time = t;
		}

		AggregateOutput ao = new AggregateOutput();
		ao.ttError = ttError;
		ao.ttTime = ttTime;
		ao.inputHash = ri.hash;
		ao.params = ri.params;
		ao.mos = mos;

		ao.outputHash = Utils.getHash(hashsb.toString());
		// as such, the hash of the final file is different

		logger.info("\nTotal error: {}\n" +
				"Single-threaded time taken: {}\n" +
				"Output hash: {}", ao.ttError, ao.ttTime, ao.outputHash);

		FileWriter fw = new FileWriter(
				output + "-" + ri.hash + " (" + ao.outputHash + ").json");
		fw.write(gson.toJson(ao));
		fw.close();
	}

	public static void outputOne(MoleculeOutput mo, String output)
			throws IOException {
		GsonBuilder builder = new GsonBuilder();
		builder.serializeNulls();
		builder.serializeSpecialFloatingPointValues();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		String jsoned = gson.toJson(mo);
		FileWriter fw = new FileWriter(output + ".json", true);
		fw.write(jsoned + ",\n");
		fw.close();
	}

	public static void main(String[] args) {
		MoleculeOutput[] ranMolecules =
				OutputHandler.importMoleculeOutputs(
						"outputs/run-0000-output-1547508A");
		long time = 0;
		long max = 0;
		double ttError = 0;
		for (MoleculeOutput mo : ranMolecules) {
			time += mo.time;
			if (mo.time > max) max = mo.time;
			ttError += mo.totalError;
			System.out.println(
					mo.rawMolecule.index + "mo.totalError = " + mo.totalError);
		}
		System.out.println(time / 1e3 / 60 / 60);
		System.out.println("max = " + max / 1e3 / 60);
		System.out.println("ttError = " + ttError);
	}
}
