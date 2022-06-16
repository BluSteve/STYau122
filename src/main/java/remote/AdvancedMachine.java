package remote;

import com.google.gson.JsonSyntaxException;
import host.Machine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import runcycle.IMoleculeResult;
import runcycle.structs.InputInfo;
import runcycle.structs.RunInput;
import runcycle.structs.RunnableMolecule;
import runcycle.structs.Serializer;
import tools.Byter;

import java.net.Socket;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import static remote.Utils.fromJsonBytes;
import static runcycle.structs.Serializer.gson;

public class AdvancedMachine extends Machine implements Comparable<AdvancedMachine> {
	public final Logger logger;
	public double power;

	public AdvancedMachine(Socket socket) {
		super(socket);
		logger = LogManager.getLogger(name);
	}

	public void init() {
		sendJar("build/libs/MNDOParam.jar", "remote.AdditionalMethods");
		power = getPower();
	}

	private void write(String methodName, Object object) {
		super.write(methodName, gson.toJson(object).getBytes(StandardCharsets.UTF_8));
	}

	private double getPower() {
		write("benchmark");
		return Byter.toDouble(read());
	}

	public void updatePower() {
		write("updatePower", Byter.toBytes(power));
	}

	public RunInput buildMolecules(String pnFile, String pFile, String mFile) {
		write("buildMolecules", new MoleculeInputFiles(pnFile, pFile, mFile));

		return fromJsonBytes(read(), RunInput.class);
	}

	public String getLogs() {
		write("getLogs");

		return new String(read());
	}

	private long[] getMoleculesHash() {
		write("getMoleculesHash");

		return fromJsonBytes(read(), long[].class);
	}

	public SubsetResult runMolecules(RunnableMolecule[] rms, InputInfo info, String inputHash) {
		long[] cacheHashes = getMoleculesHash();

		List<RunnableMolecule> toSend = new ArrayList<>();
		List<Long> toUseCached = new ArrayList<>();
		for (RunnableMolecule rm : rms) {
			long hash = Serializer.getLongHash(rm);
			boolean present = false;
			for (long cacheHash : cacheHashes) {
				if (hash == cacheHash) {
					present = true;
					break;
				}
			}
			if (present) toUseCached.add(hash);
			else toSend.add(rm); // add rms which are not in cache

		}
		RunnableMolecule[] toSendArr = toSend.toArray(new RunnableMolecule[0]);
		long[] toUseCachedArr = toUseCached.stream().mapToLong(i -> i).toArray();

		System.out.println("toSendArr.length = " + toSendArr.length);
		System.out.println("toUseCachedArr.length = " + toUseCachedArr.length);

		logger.info("Uploading {} out of {} molecules.", toSendArr.length, rms.length);
		write("runMolecules", new MoleculesSubset(toSendArr, toUseCachedArr, info, inputHash));


		List<IMoleculeResult> hithertoResults = new LinkedList<>();

		SubsetResult sr;
		while (true) {
			byte[] bytes = read();
			try {
				sr = fromJsonBytes(bytes, SubsetResult.class);
				break;
			} catch (JsonSyntaxException e) {
				IMoleculeResult[] results = fromJsonBytes(bytes, IMoleculeResult[].class);
				logger.info("Downloaded {} molecules.", results.length);
				hithertoResults.addAll(List.of(results));
			}
		}

		hithertoResults.addAll(List.of(sr.results));
		sr.results = hithertoResults.toArray(new IMoleculeResult[0]);

		return sr;
	}

	@Override
	public int compareTo(AdvancedMachine o) {
		return Double.compare(o.power, power);
	}

	@Override
	public String toString() {
		return "AdvancedMachine{" +
				"name='" + name + '\'' +
				", power=" + power +
				'}';
	}

	public static class SubsetResult {
		public long timeTaken;
		public IMoleculeResult[] results;

		public SubsetResult() {
		}

		public SubsetResult(long timeTaken, IMoleculeResult[] results) {
			this.timeTaken = timeTaken;
			this.results = results;
		}
	}

	static class MoleculeInputFiles {
		public String pnFile, pFile, mFile;

		public MoleculeInputFiles() {
		}

		public MoleculeInputFiles(String pnFile, String pFile, String mFile) {
			this.pnFile = pnFile;
			this.pFile = pFile;
			this.mFile = mFile;
		}
	}

	class MoleculesSubset {
		public RunnableMolecule[] rms;
		public long[] cachedHashes;
		public InputInfo info;
		public String inputHash, ip;

		public MoleculesSubset() {
		}

		public MoleculesSubset(RunnableMolecule[] rms, long[] cachedHashes, InputInfo info, String inputHash) {
			this.rms = rms;
			this.cachedHashes = cachedHashes;
			this.info = info;
			this.inputHash = inputHash;
			this.ip = AdvancedMachine.this.ip;
		}
	}
}
