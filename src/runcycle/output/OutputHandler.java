package runcycle.output;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import runcycle.MoleculeRun;
import runcycle.input.RawInput;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class OutputHandler {
	public static MoleculeOutput toMoleculeOutput(MoleculeRun result) {
		MoleculeOutput mo = new MoleculeOutput();
		mo.rawMolecule = result.getRm();
		mo.time = result.getTime();
		if (result.getH() != null) mo.hessian = result.getH().getHessianRaw();

		mo.hf = result.getS().hf;
		mo.dipole = result.getS().dipole;
		mo.ie = -result.getS().homo;
		mo.geomGradient = result.getE().getGeomGradient();
		mo.totalError = result.getE().getTotalError();

		ParamGradientOutput pgo = new ParamGradientOutput();
		pgo.hf = result.getG().getHFDerivs();
		pgo.dipole = result.getG().getDipoleDerivs();
		pgo.ie = result.getG().getIEDerivs();
		pgo.geom = result.getG().getGeomDerivs();
		pgo.total = result.getG().getTotalGradients();

		mo.gradient = pgo;
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
		builder.setPrettyPrinting();
		Gson gson = builder.create();
		try {
			String jsoned = gson.toJson(mo);
			FileWriter fw = new FileWriter(output + ".json", true);
			fw.write(jsoned + ",\n");
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		MoleculeOutput[] ranMolecules = OutputHandler.importMoleculeOutputs("dynamic-output");
	}
}
