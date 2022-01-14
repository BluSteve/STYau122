package runcycle.structs;

import com.google.gson.*;
import org.apache.commons.lang3.StringUtils;
import runcycle.IMoleculeResult;
import tools.Utils;

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
		InputInfo info = ro.nextInput.info;

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
			json.add("paramsMap", gson.toJsonTree(src.getParamsMap()));

			return json;
		};

		JsonDeserializer<InputInfo> iids = (json, typeOfT, context) -> {
			JsonObject object = json.getAsJsonObject();

			return new InputInfo(
					gson.fromJson(object.get("atomTypes"), int[].class),
					gson.fromJson(object.get("neededParams"), int[][].class),
					gson.fromJson(object.get("paramsMap"), double[][].class)
			);
		};

		JsonSerializer<RunnableMolecule> rms = (src, typeOfSrc, context) -> {
			JsonObject json = gson.toJsonTree(src).getAsJsonObject();

			json.remove("atoms");
			json.remove("expGeom");

			Atom[] bohratoms = new Atom[src.atoms.length];
			for (int i = 0; i < bohratoms.length; i++) {
				bohratoms[i] = new Atom(src.atoms[i].Z, Utils.debohr(src.atoms[i].coords));
			}

			json.add("atoms", gson.toJsonTree(bohratoms));

			if (src.expGeom != null) {
				Atom[] bohrexp = new Atom[src.expGeom.length];
				for (int i = 0; i < bohratoms.length; i++) {
					bohrexp[i] = new Atom(src.expGeom[i].Z, Utils.debohr(src.expGeom[i].coords));
				}
				json.add("expGeom", gson.toJsonTree(bohrexp));
			}

			return json;
		};

		JsonDeserializer<RunnableMolecule> rmds = (json, typeOfT, context) -> {
			JsonObject object = json.getAsJsonObject();

			RunnableMolecule.RMBuilder builder = new RunnableMolecule.RMBuilder();

			builder.index = object.get("index").getAsInt();
			builder.name = object.get("name").getAsString();
			builder.restricted = object.get("restricted").getAsBoolean();
			builder.useEdiis = object.get("useEdiis").getAsBoolean();
			builder.charge = object.get("charge").getAsInt();
			builder.mult = object.get("mult").getAsInt();
			builder.datum = gson.fromJson(object.get("datum"), double[].class);
			builder.atoms = gson.fromJson(object.get("atoms"), Atom[].class);
			builder.expGeom = gson.fromJson(object.get("expGeom"), Atom[].class);
			builder.densityMatrices = gson.fromJson(object.get("densityMatrices"), double[][].class);
			builder.densityMatricesExp =  gson.fromJson(object.get("densityMatricesExp"), double[][].class);

			for (int i = 0; i < builder.atoms.length; i++) {
				builder.atoms[i] = new Atom(builder.atoms[i].Z, Utils.bohr(builder.atoms[i].coords));
				if (builder.expGeom != null)
					builder.expGeom[i] = new Atom(builder.expGeom[i].Z, Utils.bohr(builder.expGeom[i].coords));
			}

			int[] mats = gson.fromJson(object.get("mats"), int[].class);
			int[][] mnps = gson.fromJson(object.get("mnps"), int[][].class);

			return builder.build(mats, mnps, null);
		};

		GsonBuilder builder = new GsonBuilder();
		builder.registerTypeAdapter(InputInfo.class, iis);
		builder.registerTypeAdapter(InputInfo.class, iids);
		builder.registerTypeAdapter(RunnableMolecule.class, rms);
		builder.registerTypeAdapter(RunnableMolecule.class, rmds);

		Gson finalgson = builder.create();

		JsonSerializer<IMoleculeResult> imrs =
				(src, typeOfSrc, context) -> finalgson.toJsonTree(MoleculeOutput.from(src));

		JsonDeserializer<IMoleculeResult> imrds =
				(json, typeOfT, context) -> finalgson.fromJson(json, MoleculeOutput.class);

		builder.registerTypeAdapter(IMoleculeResult.class, imrs);
		builder.registerTypeAdapter(IMoleculeResult.class, imrds);
//		builder.setPrettyPrinting();

		return builder.create();
	}

	private static class MoleculeOutput implements IMoleculeResult {
		public RunnableMolecule updatedRm;
		public long time;
		public double hf, dipole, ie, geomGradMag, totalError;
		public double[][] hfpg, dipolepg, iepg, geompg, totalpg, hessian;

		public MoleculeOutput() {
		}

		private MoleculeOutput(RunnableMolecule updatedRm, long time, double hf, double dipole, double ie,
							   double geomGradMag, double totalError, double[][] hfpg, double[][] dipolepg,
							   double[][] iepg, double[][] geompg, double[][] totalpg, double[][] hessian) {
			this.updatedRm = updatedRm;
			this.time = time;
			this.hf = hf;
			this.dipole = dipole;
			this.ie = ie;
			this.geomGradMag = geomGradMag;
			this.totalError = totalError;
			this.hfpg = hfpg;
			this.dipolepg = dipolepg;
			this.iepg = iepg;
			this.geompg = geompg;
			this.totalpg = totalpg;
			this.hessian = hessian;
		}

		public static MoleculeOutput from(IMoleculeResult result) {
			return new MoleculeOutput(result.getUpdatedRm(), result.getTime(), result.getHf(),
					result.getDipole(), result.getIE(), result.getGeomGradMag(), result.getTotalError(),
					result.getHfDerivs(), result.getDipoleDerivs(), result.getIEDerivs(), result.getGeomDerivs(),
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
		public double getHf() {
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
		public double getGeomGradMag() {
			return geomGradMag;
		}

		@Override
		public double getTotalError() {
			return totalError;
		}

		@Override
		public double[][] getHfDerivs() {
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
