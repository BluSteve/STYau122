package scf;

import java.io.Serializable;
import java.util.Arrays;

public class AtomProperties implements Serializable { // Only 119 of these
	private static final double HEATCONV = 4.3363E-2;
	private int index, Z, Q;
	private String name;
	private double mass, heat;
	private OrbitalProperties[] orbitals;

	public int getPeriod() {
		int period = 1;
		if (Z > 85) period = 7;
		else if (Z > 53) period = 6;
		else if (Z > 35) period = 5;
		else if (Z > 17) period = 4;
		else if (Z > 9) period = 3;
		else if (Z > 1) period = 2;
		return period;
	}

	public int getZ() {
		return Z;
	}

	public int getQ() {
		return Q;
	}

	public double getHeat() {
		return heat * HEATCONV;
	}

	public double getMass() {
		return mass;
	}

	public String getName() {
		return name;
	}

	public OrbitalProperties[] getOrbitals() {
		return orbitals;
	}

	public void setOrbitals(OrbitalProperties[] orbitals) {
		this.orbitals = orbitals;
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}

	@Override
	public String toString() {
		return "AtomProperties{" +
				"index=" + index +
				", Z=" + Z +
				", Q=" + Q +
				", name='" + name + '\'' +
				", mass=" + mass +
				", heat=" + heat +
				", orbitals=" + Arrays.toString(orbitals) +
				'}';
	}
}
