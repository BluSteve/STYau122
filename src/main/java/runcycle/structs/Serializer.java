package runcycle.structs;

import com.google.gson.*;
import org.apache.commons.lang3.StringUtils;
import runcycle.IMoleculeResult;

import java.nio.ByteBuffer;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

public class Serializer {
	public static final Gson gson = getGson();

	public static String getHash(String str) {
		try {
			MessageDigest digest = MessageDigest.getInstance("SHA-1");
			byte[] b = digest.digest(str.getBytes(StandardCharsets.UTF_8));
			long v = ByteBuffer.wrap(b).getLong();

			// gets last 51 bits of hash, 36^10 is 51.29 bits of hash
			// ensures low collision probability up to 10k runs
			return StringUtils.leftPad(Long.toUnsignedString(v & 0x7FFFFFFFFFFFFL, 36).toUpperCase(), 10, "0");
		} catch (NoSuchAlgorithmException e) {
			throw new RuntimeException(e);
		}
	}

	public static String getHash(RunOutput ro) {
		MoleculeOutput[] mos = new MoleculeOutput[ro.results.length];
		for (int i = 0; i < mos.length; i++) {
			mos[i] = MoleculeOutput.from(ro.results[i]);
			mos[i].time = 0;
		}
		InputInfo info = ro.nextRunInfo;

		return getHash(gson.toJson(info) + gson.toJson(mos));
	}

	public static String getHash(RunInput ri) {
		return getHash(gson.toJson(ri));
	}

	public static String getHash(Object o) {
		return getHash(gson.toJson(o));
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
				(src, typeOfSrc, context) -> gson.toJsonTree(MoleculeOutput.from(src));

		JsonDeserializer<IMoleculeResult> imrds =
				(json, typeOfT, context) -> gson.fromJson(json, MoleculeOutput.class);

		GsonBuilder builder = new GsonBuilder();
		builder.registerTypeAdapter(InputInfo.class, iis);
		builder.registerTypeAdapter(InputInfo.class, iids);
		builder.registerTypeAdapter(RunnableMolecule.class, rms);
		builder.registerTypeAdapter(RunnableMolecule.class, rmds);
		builder.registerTypeAdapter(IMoleculeResult.class, imrs);
		builder.registerTypeAdapter(IMoleculeResult.class, imrds);
//		builder.setPrettyPrinting();

		return builder.create();
	}

	private static class MoleculeOutput implements IMoleculeResult {
		public RunnableMolecule updatedRm;
		public long time;
		public double hf, dipole, ie, geomGradient, totalError;
		public double[][] hfpg, dipolepg, iepg, geompg, totalpg, hessian;

		public MoleculeOutput() {
		}

		private MoleculeOutput(RunnableMolecule updatedRm, long time, double hf, double dipole, double ie,
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

		public static MoleculeOutput from(IMoleculeResult result) {
			return new MoleculeOutput(result.getUpdatedRm(), result.getTime(), result.getHF(),
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
