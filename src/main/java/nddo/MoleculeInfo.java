package nddo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class MoleculeInfo { // todo check if moleculeinfo is valid before using it
	public int index; // optional, for debugging
	public String name;
	public boolean restricted;
	public int charge, mult;
	public double[] datum;
	public int[] atomicNumbers;
	public int nElectrons, nOrbitals, nIntegrals, nCoulombInts, nExchangeInts;
	public int[][] orbsOfAtom, missingOfAtom;
	public int[] atomOfOrb;
	public int[] mats; // molecule atom types
	public int[][] mnps; // molecule needed params

	private transient String debugName;
	private transient Logger logger;

	public String debugName() {
		if (debugName == null)
			debugName = String.format("%03d-%s", index, name);

		return debugName;
	}

	public Logger getLogger() {
		if (logger == null) {
			logger = LogManager.getLogger(debugName());
		}

		return logger;
	}
}
