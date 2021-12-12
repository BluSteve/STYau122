package scf;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;

public class OrbitalProperties implements Serializable {
	private static final String[][] ORBITALS_PER_SHELL =
			new String[][]{new String[0], new String[]{"s"},
					new String[]{"s", "p"},
					new String[]{"s", "p", "d"},
					new String[]{"s", "p", "d", "f"},
					new String[]{"s", "p", "d", "f"},
					new String[]{"s", "p", "d", "f"},
					new String[]{"s", "p", "d", "f"}};
	private final String type;
	private final int shell, L;
	private final int[] config;

	private OrbitalProperties(String type, int shell, int[] config) {
		this.type = type;
		this.shell = shell;
		this.config = config;
		this.L = config[0] + config[1] + config[2];
	}

	static OrbitalProperties[] generateOrbitals(int shell) {
		ArrayList<OrbitalProperties> orbitals = new ArrayList<>(25);

		for (String x : ORBITALS_PER_SHELL[shell]) {
			switch (x) {
				case "s":
					orbitals.add(new OrbitalProperties(x, shell,
							new int[]{0, 0, 0}));
					break;
				case "p":
					orbitals.add(new OrbitalProperties(x, shell,
							new int[]{1, 0, 0}));
					orbitals.add(new OrbitalProperties(x, shell,
							new int[]{0, 1, 0}));
					orbitals.add(new OrbitalProperties(x, shell,
							new int[]{0, 0, 1}));
					break;
			}
		}

		return orbitals.toArray(new OrbitalProperties[0]);
	}

	public String getType() {
		return type;
	}

	public int getShell() {
		return shell;
	}

	public int[] getConfig() {
		return config;
	}

	public int getL() {
		return L;
	}

	@Override
	public String toString() {
		return "OrbitalProperties{" +
				"type='" + type + '\'' +
				", shell=" + shell +
				", L=" + L +
				", config=" + Arrays.toString(config) +
				'}';
	}
}
