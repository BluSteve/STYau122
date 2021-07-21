package runcycle.output;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import runcycle.MoleculeResult;
import runcycle.input.RawInput;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class OutputHandler {
	public static MoleculeOutput toMoleculeOutput(MoleculeResult result) {
		MoleculeOutput mo = new MoleculeOutput();
		mo.rawMolecule = result.getRm();
		mo.time = result.getTime();
		mo.isExpAvail = result.isExpAvail();

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

		if (result.getHessian() != null) mo.hessian = result.getHessian();

		return mo;
	}

	public static MoleculeOutput[] importMoleculeOutputs(String inputPath) {
		MoleculeOutput[] mos = null;
		try {
			mos = new Gson().fromJson(new FileReader(inputPath + ".json"),
					MoleculeOutput[].class);
		} catch (Exception ignored) {
		}
		return mos;
	}

	public static void output(RawInput ri, MoleculeOutput[] mos,
							  String output) {
		GsonBuilder builder = new GsonBuilder();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		try {
			String jsoned = gson.toJson(mos);
			double ttError = 0;
			for (MoleculeOutput mo : mos) {
				ttError += mo.totalError;
			}
			System.err.println("\nTotal Error: " + ttError);

			// uses the hash of input to distinguish outputs. if the program is
			// working correctly inputs should map to outputs one-to-one.
			FileWriter fw = new FileWriter(output + "-" + ri.hash + ".json");
			fw.write(jsoned);
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void outputOne(MoleculeOutput mo, String output) {
		GsonBuilder builder = new GsonBuilder();
		builder.serializeNulls();
		builder.serializeSpecialFloatingPointValues();
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		try {
			String jsoned = gson.toJson(mo);
			FileWriter fw = new FileWriter(output + ".json", true);
			fw.write(jsoned + ",\n");
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			e.printStackTrace();
			System.err.println(mo.toString());
		}
	}

	public static void main(String[] args) {
		MoleculeOutput[] ranMolecules =
				OutputHandler.importMoleculeOutputs("dynamic-output");
		long time = 0;
		long max = 0;
		for (MoleculeOutput mo : ranMolecules) {
			time += mo.time;
			if (mo.time > max) max = mo.time;
		}
		System.out.println(time / 1e3 / 60 / 60);
		System.out.println("max = " + max / 1e3 / 60);
	}
}
