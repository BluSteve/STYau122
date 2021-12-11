package runcycle.input;

import nddo.MoleculeInfo;

import java.util.Arrays;

public class RawMolecule extends MoleculeInfo {
	public RawAtom[] atoms, expGeom;

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
				'}';
	}
}
