package runcycle.input;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Arrays;

public class RawMolecule {
	public int index;
	public String name;
	public boolean restricted;
	public int charge, mult;
	public double[] datum;
	public int[] atomicNumbers;
	public int nElectrons, nOrbitals, nIntegrals, nCoulombInts, nExchangeInts;
	public int[] mats; // molecule atom types
	public int[][] mnps; // molecule needed params
	public RawAtom[] atoms, expGeom;
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

	@Override
	public String toString() {
		return "RawMolecule{" +
				"index=" + index +
				", name='" + name + '\'' +
				", restricted=" + restricted +
				", charge=" + charge +
				", mult=" + mult +
				", datum=" + Arrays.toString(datum) +
				", atomicNumbers=" + Arrays.toString(atomicNumbers) +
				", nElectrons=" + nElectrons +
				", nOrbitals=" + nOrbitals +
				", nIntegrals=" + nIntegrals +
				", nCoulombInts=" + nCoulombInts +
				", nExchangeInts=" + nExchangeInts +
				", mats=" + Arrays.toString(mats) +
				", mnps=" + Arrays.toString(mnps) +
				", atoms=" + Arrays.toString(atoms) +
				", expGeom=" + Arrays.toString(expGeom) +
				", debugName='" + debugName + '\'' +
				'}';
	}
}
