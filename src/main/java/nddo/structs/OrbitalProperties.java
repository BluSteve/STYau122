package nddo.structs;

import java.util.ArrayList;
import java.util.List;

public class OrbitalProperties {
	private static final String[][] ORBITALS_PER_SHELL =
			new String[][]{new String[0], new String[]{"s"},
					new String[]{"s", "p"},
					new String[]{"s", "p", "d"},
					new String[]{"s", "p", "d", "f"},
					new String[]{"s", "p", "d", "f"},
					new String[]{"s", "p", "d", "f"},
					new String[]{"s", "p", "d", "f"}};
	private static final List<OrbitalProperties> cache = new ArrayList<>();

	public final String type;
	public final int shell, L, i, j, k;

	private OrbitalProperties(String type, int shell, int i, int j, int k) {
		this.type = type;
		this.shell = shell;
		this.i = i;
		this.j = j;
		this.k = k;
		this.L = i + j + k;
	}

	/**
	 * Implements the Flyweight pattern. This method ensures a finite number of OrbitalProperties objects.
	 *
	 * @return Potentially cached OrbitalProperties object.
	 */
	private static OrbitalProperties of(String type, int shell, int i, int j, int k) {
		for (OrbitalProperties op : cache) {
			if (op.type.equals(type) && op.shell == shell && op.i == i && op.j == j && op.k == k) {
				return op;
			}
		}

		OrbitalProperties newop = new OrbitalProperties(type, shell, i, j, k);
		cache.add(newop);

		return newop;
	}

	static OrbitalProperties[] generateOrbitals(int shell) {
		ArrayList<OrbitalProperties> orbitals = new ArrayList<>(25);

		for (String x : ORBITALS_PER_SHELL[shell]) {
			switch (x) {
				case "s":
					orbitals.add(OrbitalProperties.of(x, shell, 0, 0, 0));
					break;
				case "p":
					orbitals.add(OrbitalProperties.of(x, shell, 1, 0, 0));
					orbitals.add(OrbitalProperties.of(x, shell, 0, 1, 0));
					orbitals.add(OrbitalProperties.of(x, shell, 0, 0, 1));
					break;
				default:
					throw new IllegalArgumentException("Unidentified orbital type!");
			}
		}

		return orbitals.toArray(new OrbitalProperties[0]);
	}

	@Override
	public String toString() {
		return "OrbitalProperties{" +
				"type='" + type + '\'' +
				", shell=" + shell +
				", L=" + L +
				", i=" + i +
				", j=" + j +
				", k=" + k +
				'}';
	}
}
