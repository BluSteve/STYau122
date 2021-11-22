package runcycle.input;

import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.appender.FileAppender;
import org.apache.logging.log4j.core.layout.PatternLayout;
import scf.AtomHandler;

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

	public static MNDOAtom[] toMNDOAtoms(RawAtom[] atoms,
										 MNDOParams[] mndoParams) {
		MNDOAtom[] res = new MNDOAtom[atoms.length];
		for (int i = 0; i < atoms.length; i++) {
			res[i] = atoms[i].toMNDOAtom(
					mndoParams[AtomHandler.atoms[atoms[i].Z].getIndex()]);
		}
		return res;
	}

	public String debugName() {
		if (debugName == null)
			debugName = String.format("%03d-%s", index, name);

		return debugName;
	}

	public Logger getLogger() {
		if (logger == null) {
			logger = LogManager.getLogger(debugName());
			org.apache.logging.log4j.core.Logger l =
					(org.apache.logging.log4j.core.Logger) logger;
			LoggerContext lc = l.getContext();
			FileAppender fa = FileAppender.newBuilder()
					.setName(debugName)
					.withAppend(false)
					.withFileName("logs/" + debugName + ".log")
					.setLayout(
							PatternLayout.newBuilder().withPattern(
											"%d{ISO8601} %08r %-5level " +
													"%logger{36} - %msg%n")
									.build())
					.setConfiguration(lc.getConfiguration()).build();
			fa.start();
			lc.getConfiguration().addAppender(fa);
			((org.apache.logging.log4j.core.Logger) logger).addAppender(
					lc.getConfiguration().getAppender(fa.getName()));
			lc.updateLoggers();
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
