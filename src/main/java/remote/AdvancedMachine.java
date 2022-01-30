package remote;

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

import static remote.Utils.fromJsonBytes;
import static runcycle.structs.Serializer.gson;

public class AdvancedMachine extends Machine implements Comparable<AdvancedMachine> {
	public final Logger logger;
	public double power;

	public AdvancedMachine(Socket socket) {
		super(socket);
		logger = LogManager.getLogger(name);
		sendJar("build/libs/MNDOParam.jar");
		power = getPower();
	}

	private double getPower() {
		write("benchmark");
		return Byter.toDouble(read());
	}

	private void write(String methodName, Object object) {
		super.write(methodName, gson.toJson(object).getBytes(StandardCharsets.UTF_8));
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

	private String getMoleculesHash() {
		write("getMoleculesHash");

		return new String(read());
	}

	public IMoleculeResult[] downloadMolecules() {
		write("downloadMolecules");

		return fromJsonBytes(read(), IMoleculeResult[].class);
	}

	public SubsetResult runMolecules(RunnableMolecule[] rms, InputInfo info, String inputHash) {
		String hash = Serializer.getHash(rms);
		String cacheHash = getMoleculesHash();

		if (hash.equals(cacheHash)) {
			logger.info("Using cached runnable molecules.");
			write("runMolecules", new MoleculesSubset(info, inputHash));
		}
		else {
			logger.info("Uploading molecules...");
			write("runMolecules", new MoleculesSubset(rms, info, inputHash));
		}

		return fromJsonBytes(read(), SubsetResult.class);
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

	public static class MoleculeInputFiles {
		public String pnFile, pFile, mFile;

		public MoleculeInputFiles() {
		}

		public MoleculeInputFiles(String pnFile, String pFile, String mFile) {
			this.pnFile = pnFile;
			this.pFile = pFile;
			this.mFile = mFile;
		}
	}

	public class MoleculesSubset {
		public RunnableMolecule[] rms;
		public InputInfo info;
		public String inputHash, ip;

		public MoleculesSubset() {
		}

		public MoleculesSubset(RunnableMolecule[] rms, InputInfo info, String inputHash) {
			this.rms = rms;
			this.info = info;
			this.inputHash = inputHash;
			this.ip = AdvancedMachine.this.ip;
		}

		public MoleculesSubset(InputInfo info, String inputHash) {
			this(null, info, inputHash);
		}
	}
}
