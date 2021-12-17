package frontend.json;

import com.google.gson.*;
import runcycle.IMoleculeResult;
import runcycle.structs.Atom;
import runcycle.structs.InputInfo;
import runcycle.structs.RunnableMolecule;

import java.io.FileWriter;
import java.io.IOException;

public class Serializer {
	private static final Gson gson = getGson();

	public static void write(Object o, String first, String... filenames) throws IOException {
		FileWriter fw = new FileWriter(first + ".json");
		gson.toJson(o, fw);
		fw.close();

		for (String filename : filenames) {
			fw = new FileWriter(filename + ".json");
			gson.toJson(o, fw);
			fw.close();
		}
	}

	private static Gson getGson() {
		Gson gson = new Gson();
		JsonSerializer<InputInfo> iis = (src, typeOfSrc, context) -> {
			JsonObject json = new JsonObject();

			json.add("atomTypes", gson.toJsonTree(src.atomTypes));
			json.add("neededParams", gson.toJsonTree(src.neededParams));
			json.add("params", gson.toJsonTree(src.getParams()));

			return json;
		};

		JsonDeserializer<InputInfo> iids = (json, typeOfT, context) -> {
			JsonObject object = json.getAsJsonObject();

			return new InputInfo(
					gson.fromJson(object.get("atomTypes"), int[].class),
					gson.fromJson(object.get("neededParams"), int[][].class),
					gson.fromJson(object.get("params"), double[][].class)
			);
		};

		JsonSerializer<RunnableMolecule> rms = (src, typeOfSrc, context) -> {
			JsonObject json = new JsonObject();

			json.addProperty("index", src.index);
			json.addProperty("name", src.name);
			json.addProperty("restricted", src.restricted);
			json.addProperty("charge", src.charge);
			json.addProperty("mult", src.mult);
			json.add("datum", gson.toJsonTree(src.datum));
			json.add("mats", gson.toJsonTree(src.mats));
			json.add("mnps", gson.toJsonTree(src.mnps));
			json.add("atoms", gson.toJsonTree(src.atoms));
			json.add("expGeom", gson.toJsonTree(src.expGeom));

			return json;
		};

		JsonDeserializer<RunnableMolecule> rmds = (json, typeOfT, context) -> {
			JsonObject object = json.getAsJsonObject();

			RunnableMolecule.RMBuilder builder = new RunnableMolecule.RMBuilder();

			builder.index = object.get("index").getAsInt();
			builder.name = object.get("name").getAsString();
			builder.restricted = object.get("restricted").getAsBoolean();
			builder.charge = object.get("charge").getAsInt();
			builder.mult = object.get("mult").getAsInt();
			builder.datum = gson.fromJson(object.get("datum"), double[].class);
			builder.atoms = gson.fromJson(object.get("atoms"), Atom[].class);
			builder.expGeom = gson.fromJson(object.get("expGeom"), Atom[].class);

			int[] mats = gson.fromJson(object.get("mats"), int[].class);
			int[][] mnps = gson.fromJson(object.get("mnps"), int[][].class);

			return builder.build(mats, mnps);
		};

		JsonSerializer<IMoleculeResult> imrs =
				(src, typeOfSrc, context) -> gson.toJsonTree(MoleculeSerial.from(src));

		JsonDeserializer<IMoleculeResult> imrds =
				(json, typeOfT, context) -> gson.fromJson(json, MoleculeSerial.class);

		GsonBuilder builder = new GsonBuilder();
		builder.registerTypeAdapter(InputInfo.class, iis);
		builder.registerTypeAdapter(InputInfo.class, iids);
		builder.registerTypeAdapter(RunnableMolecule.class, rms);
		builder.registerTypeAdapter(RunnableMolecule.class, rmds);
		builder.registerTypeAdapter(IMoleculeResult.class, imrs);
		builder.registerTypeAdapter(IMoleculeResult.class, imrds);
		builder.setPrettyPrinting();

		return builder.create();
	}

	private static class MoleculeSerial implements IMoleculeResult {
		public RunnableMolecule updatedRm;
		public long time;
		public double hf, dipole, ie, geomGradient, totalError;
		public double[][] hfpg, dipolepg, iepg, geompg, totalpg, hessian;

		public MoleculeSerial() {
		}

		private MoleculeSerial(RunnableMolecule updatedRm, long time, double hf, double dipole, double ie,
							   double geomGradient, double totalError, double[][] hfpg, double[][] dipolepg,
							   double[][] iepg, double[][] geompg, double[][] totalpg, double[][] hessian) {
			this.updatedRm = updatedRm;
			this.time = time;
			this.hf = hf;
			this.dipole = dipole;
			this.ie = ie;
			this.geomGradient = geomGradient;
			this.totalError = totalError;
			this.hfpg = hfpg;
			this.dipolepg = dipolepg;
			this.iepg = iepg;
			this.geompg = geompg;
			this.totalpg = totalpg;
			this.hessian = hessian;
		}

		public static MoleculeSerial from(IMoleculeResult result) {
			return new MoleculeSerial(result.getUpdatedRm(), result.getTime(), result.getHF(),
					result.getDipole(), result.getIE(), result.getGeomGradient(), result.getTotalError(),
					result.getHFDerivs(), result.getDipoleDerivs(), result.getIEDerivs(), result.getGeomDerivs(),
					result.getTotalGradients(), result.getHessian());
		}

		@Override
		public RunnableMolecule getUpdatedRm() {
			return updatedRm;
		}

		@Override
		public long getTime() {
			return time;
		}

		@Override
		public boolean isExpAvail() {
			return updatedRm.expGeom != null;
		}

		@Override
		public double[][] getHessian() {
			if (hessian != null) return hessian;
			else throw new IllegalStateException(
					"Hessian not found for previously ran molecule: " + updatedRm.debugName());
		}

		@Override
		public double getHF() {
			return hf;
		}

		@Override
		public double getDipole() {
			return dipole;
		}

		@Override
		public double getIE() {
			return ie;
		}

		@Override
		public double getGeomGradient() {
			return geomGradient;
		}

		@Override
		public double getTotalError() {
			return totalError;
		}

		@Override
		public double[][] getHFDerivs() {
			return hfpg;
		}

		@Override
		public double[][] getDipoleDerivs() {
			return dipolepg;
		}

		@Override
		public double[][] getIEDerivs() {
			return iepg;
		}

		@Override
		public double[][] getGeomDerivs() {
			return geompg;
		}

		@Override
		public double[][] getTotalGradients() {
			return totalpg;
		}
	}
}
