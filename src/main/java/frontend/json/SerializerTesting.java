package frontend.json;

import com.google.gson.*;
import frontend.txt.Main;
import runcycle.IMoleculeResult;
import runcycle.RunIterator;
import runcycle.structs.*;

import java.io.IOException;

public class SerializerTesting {
	public static void main(String[] args) throws IOException {
		InputInfo info = Main.getInputInfo();
		RunnableMolecule[] molecules = Main.getMolecules(info.atomTypes, info.neededParams);
		RunIterator ri = new RunIterator(new RunInput(info, molecules));
		IMoleculeResult imr = ri.next().results[0];

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

		Gson g = builder.create();
//		System.out.println(g.toJson(info));
//		InputInfo i = g.fromJson(g.toJson(info), InputInfo.class);
//		RunnableMolecule rm = g.fromJson(g.toJson(imr.getUpdatedRm()), RunnableMolecule.class);
//		IMoleculeResult imr2 = g.fromJson(g.toJson(MoleculeSerial.from(imr)), MoleculeSerial.class);
		RunOutput ro = g.fromJson(g.toJson(ri.next()), RunOutput.class);
		System.exit(0);
	}
}
