package scf;

import org.apache.logging.log4j.LogManager;

import java.util.Objects;

public class STO6G extends LCGTO {

	private final static double[] exp1 =
			new double[]{6.51095E-2, 1.58088E-1, 4.07099E-1, 1.18506, 4.23592,
					2.31030E+1};

	private final static double[] coeff1 =
			new double[]{1.30334E-1, 4.16492E-1, 3.70563E-1, 1.68538E-1,
					4.93615E-2,
					9.16360E-3};

	private final static double[] exp2 =
			new double[]{4.85690E-2, 1.05960E-1, 2.43977E-1, 6.34142E-1,
					2.04036,
					1.03087E+1};

	private final static double[] coeff2s =
			new double[]{2.40706E-1, 5.95117E-1, 2.50242E-1, -3.37854E-2,
					-4.69917E-2,
					-1.32528E-2};

	private final static double[] coeff2p =
			new double[]{1.01708E-1, 4.25860E-1, 4.18036E-1, 1.73897E-1,
					3.76794E-2,
					3.75970E-3};

	public STO6G(double zeta, AtomFixed atom, OrbitalProperties orbital) {
		super(exp(orbital.getShell(), zeta), Objects.requireNonNull(coeff(orbital.getShell(), orbital.getL())), atom,
				orbital, zeta);
	}

	public STO6G() {

	}

	private static double[] exp(int shell, double zeta) {
		switch (shell) {
			case 1:
				double[] exp1 = new double[6];
				for (int i = 0; i < 6; i++) {
					exp1[i] = STO6G.exp1[i] * zeta * zeta;
				}
				return exp1;
			case 2:
				double[] exp2 = new double[6];
				for (int i = 0; i < 6; i++) {
					exp2[i] = STO6G.exp2[i] * zeta * zeta;
				}
				return exp2;
		}

		LogManager.getLogger().error("Illegal atom");
		return null;
	}

	private static double[] coeff(int shell, int L) {
		switch (shell) {
			case 1:
				return coeff1.clone();
			case 2:
				switch (L) {
					case 0:
						return coeff2s.clone();
					case 1:
						return coeff2p.clone();
				}
		}

		LogManager.getLogger().error("Illegal atom" + shell + L);
		return null;
	}
}
