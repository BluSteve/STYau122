package mndoparam.mndo;

import scf.*;

public class MNDOAtom extends AtomFixed {
    //    protected double mp.alpha, mp.betas, mp.betap, mp.Uss, mp.Upp, mp.zetas, mp.zetap, mp.eisol, mp.gss, mp.gsp, mp.hsp, mp.gpp, mp.gp2
    protected double p0, p1, p2, D1, D2;
    private MNDO6G[] orbitals;
    private MNDOParams mp;

    // TODO make MNDOParams

    public MNDOAtom(AtomProperties atomProperties, double[] coordinates, MNDOParams mp) {
        super(atomProperties, coordinates);
        this.mp = mp;
        this.p0 = p0();
        this.D1 = D1();
        this.D2 = D2();
        this.p1 = p1();
        this.p2 = p2();
        this.orbitals = setOrbitals();
    }

    //public MNDOAtom(double[] coords, MNDOAtom a) {
    public MNDOAtom(MNDOAtom a, double[] coords) {
        this(a.atomProperties, coords.clone(), a.mp);
    }

    public MNDOAtom(MNDOAtom a) {
        this(a.atomProperties, a.coordinates.clone(), a.mp.clone());
    }

    public MNDOAtom(MNDOAtom a, MNDOParams mp) {
        this(a.atomProperties, a.coordinates.clone(), mp.clone());
    }

    public MNDO6G[] getOrbitals() {
        return this.orbitals;
    }

    private MNDO6G[] setOrbitals() {//initialises OM2-3G basis functions
        OrbitalProperties[] orbitalProperties = super.atomProperties.getOrbitals();
        MNDO6G[] mndoOrbitals = new MNDO6G[orbitalProperties.length];
        for (int x = 0; x < mndoOrbitals.length; x++) {
            switch (orbitalProperties[x].getType()) {
                case "s":
                    mndoOrbitals[x] = new MNDO6G(this, orbitalProperties[x], mp.getZetas(), mp.getBetas(), mp.getUss());
                    break;
                case "p":
                    mndoOrbitals[x] = new MNDO6G(this, orbitalProperties[x], mp.getZetap(), mp.getBetap(), mp.getUpp());
                    break;
            }
        }
        return mndoOrbitals;
    }

    public double getMass() {
        return AtomHandler.atoms[atomProperties.getZ()].getMass();
    }

    public double getHeat() {
        return AtomHandler.atoms[atomProperties.getZ()].getHeat();
    }

    public String getName() {
        return AtomHandler.atoms[atomProperties.getZ()].getName();
    }

    public MNDOParams getParams() {
        return mp.clone();
    }

    public MNDO6G s() {
        return this.orbitals[0];
    }

    public double V(MNDO6G a, MNDO6G b) {
        return -this.atomProperties.getQ() * MNDO6G.getG(a, b, this.s(), this.s());
    }

    public static double crf(MNDOAtom a, MNDOAtom b) {
        double f;
        double R = GTO.R(a.getCoordinates(), b.getCoordinates()) / 1.88973;
        if ((a.atomProperties.getZ() == 7 || a.atomProperties.getZ() == 8) && b.atomProperties.getZ() == 1) {
            f = 1 + R * Math.exp(-a.mp.getAlpha() * R) + Math.exp(-b.mp.getAlpha() * R);
        } else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && a.atomProperties.getZ() == 1) {
            f = 1 + R * Math.exp(-b.mp.getAlpha() * R) + Math.exp(-a.mp.getAlpha() * R);
        } else {
            f = 1 + Math.exp(-b.mp.getAlpha() * R) + Math.exp(-a.mp.getAlpha() * R);
        }

        return f * a.atomProperties.getQ() * b.atomProperties.getQ() * MNDO6G.getG(a.s(), a.s(), b.s(), b.s());
    }

    public static double crfDeriv(MNDOAtom a, MNDOAtom b, int tau) {
        double f;
        double fprime;
        double R = GTO.R(a.getCoordinates(), b.getCoordinates());
        if ((a.atomProperties.getZ() == 7 || a.atomProperties.getZ() == 8) && b.atomProperties.getZ() == 1) {
            f = 1 + R / 1.88973 * Math.exp(-a.mp.getAlpha() * R / 1.88973) + Math.exp(-b.mp.getAlpha() * R / 1.88973);

            fprime = (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * 1.88973) * Math.exp(-a.mp.getAlpha() * R / 1.88973)
                    - a.mp.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (1.88973 * 1.88973) * Math.exp(-a.mp.getAlpha() * R / 1.88973)
                    - b.mp.getAlpha() / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-b.mp.getAlpha() * R / 1.88973);
        } else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && a.atomProperties.getZ() == 1) {
            f = 1 + R / 1.88973 * Math.exp(-b.mp.getAlpha() * R / 1.88973) + Math.exp(-a.mp.getAlpha() * R / 1.88973);

            fprime = (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * 1.88973) * Math.exp(-b.mp.getAlpha() * R / 1.88973)
                    - b.mp.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (1.88973 * 1.88973) * Math.exp(-b.mp.getAlpha() * R / 1.88973)
                    - a.mp.getAlpha() / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-a.mp.getAlpha() * R / 1.88973);
        } else {
            f = 1 + Math.exp(-b.mp.getAlpha() * R / 1.88973) + Math.exp(-a.mp.getAlpha() * R / 1.88973);

            fprime = -b.mp.getAlpha() / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-b.mp.getAlpha() * R / 1.88973)
                    - a.mp.getAlpha() / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-a.mp.getAlpha() * R / 1.88973);
        }

        return fprime * a.atomProperties.getQ() * b.atomProperties.getQ() * MNDO6G.getG(a.s(), a.s(), b.s(), b.s()) + f * a.atomProperties.getQ() * b.atomProperties.getQ() * MNDODerivative.getGderiv(a.s(), a.s(), b.s(), b.s(), tau);
    }

    public double getEisol() {
        return mp.getEisol();
    }

    private double p0() {
        return 27.2114 / (2 * mp.getGss());
    }

    private double D1() {
        return (2 * atomProperties.getPeriod() + 1) / Math.sqrt(3) * Math.pow(4 * mp.getZetas() * mp.getZetap(), atomProperties.getPeriod() + 0.5) / Math.pow(mp.getZetas() + mp.getZetap(), 2 * atomProperties.getPeriod() + 2);
    }

    private double D2() {
        return 1 / mp.getZetap() * Math.sqrt((2 * atomProperties.getPeriod() + 1) * (2 * atomProperties.getPeriod() + 2) / 20.0);
    }

    private double p1() {
        double guess = 0;

        double newguess = 0.5 * Math.pow(D1 * D1 * 27.2114 / (mp.getHsp()), 1.0 / 3);

        while (Math.abs(guess - newguess) > 1E-12) {

            guess = newguess;
            double f = 1 / guess - 1 / Math.sqrt(guess * guess + D1 * D1) - 4 * mp.getHsp() / 27.2114;
            double fprime = -1 / (guess * guess) + guess / Math.pow(guess * guess + D1 * D1, 1.5);

            newguess = guess - f / fprime;
        }
        return newguess;
    }

    private double p2() {
        double guess = 0;
        double newguess = 0.5;

        while (Math.abs(guess - newguess) > 1E-12) {
            guess = newguess;
            double f = 1 / guess + 1 / Math.sqrt(guess * guess + 2 * D2 * D2) - 2 / Math.sqrt(guess * guess + D2 * D2) - 4 * (mp.getGpp() - mp.getGp2()) / 27.2114;
            double fprime = -1 / (guess * guess) - guess / Math.pow(guess * guess + 2 * D2 * D2, 1.5) + 2 * guess / Math.pow(guess * guess + D2 * D2, 1.5);

            newguess = guess - f / fprime;
        }
        return newguess;
    }


}
