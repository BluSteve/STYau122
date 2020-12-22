package nddoparam.mndo;

import nddoparam.NDDO6G;
import nddoparam.NDDOAtom;
import scf.*;

public class MNDOAtom extends NDDOAtom {
    private MNDOParams mp;

    public MNDOAtom(AtomProperties atomProperties, double[] coordinates, MNDOParams mp) {
        super(atomProperties, coordinates, mp.nddoParams);
        this.mp = mp;
    }

    // assign new coords
    public MNDOAtom(MNDOAtom a, double[] coords) {
        this(a.atomProperties, coords.clone(), a.mp);
    }

    public MNDOAtom(MNDOAtom a) { // copy constructor
        this(a.atomProperties, a.coordinates.clone(), a.mp.clone());
    }

    // assign new params
    public MNDOAtom(MNDOAtom a, MNDOParams mp) {
        this(a.atomProperties, a.coordinates.clone(), mp.clone());
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

        return f * a.atomProperties.getQ() * b.atomProperties.getQ() * NDDO6G.getG(a.s(), a.s(), b.s(), b.s());
    }

    public static double crfDeriv(MNDOAtom a, MNDOAtom b, int tau) {
        double f;
        double fprime;
        double R = GTO.R(a.getCoordinates(), b.getCoordinates());
        if ((a.atomProperties.getZ() == 7 || a.atomProperties.getZ() == 8) && b.atomProperties.getZ() == 1) { // TODO duplicate code
            f = getF(a, b, R);
            fprime = getFPrime(a, b, R, tau);
        } else if ((b.atomProperties.getZ() == 7 || b.atomProperties.getZ() == 8) && a.atomProperties.getZ() == 1) {
            f = getF(b, a, R);
            fprime = -getFPrime(b, a, R, tau);
        } else {
            f = 1 + Math.exp(-b.mp.getAlpha() * R / 1.88973) + Math.exp(-a.mp.getAlpha() * R / 1.88973);

            fprime = -b.mp.getAlpha() / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-b.mp.getAlpha() * R / 1.88973)
                    - a.mp.getAlpha() / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-a.mp.getAlpha() * R / 1.88973);
        }

        return fprime * a.atomProperties.getQ() * b.atomProperties.getQ() * NDDO6G.getG(a.s(), a.s(), b.s(), b.s()) + f * a.atomProperties.getQ() * b.atomProperties.getQ() * MNDODerivative.getGderiv(a.s(), a.s(), b.s(), b.s(), tau);
    }

    private static double getF(MNDOAtom a, MNDOAtom b, double R) {
        return 1 + R / 1.88973 * Math.exp(-a.mp.getAlpha() * R / 1.88973) + Math.exp(-b.mp.getAlpha() * R / 1.88973);
    }

    private static double getFPrime(MNDOAtom a, MNDOAtom b, double R, int tau) {
        return (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (R * 1.88973) * Math.exp(-a.mp.getAlpha() * R / 1.88973)
                - a.mp.getAlpha() * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / (1.88973 * 1.88973) * Math.exp(-a.mp.getAlpha() * R / 1.88973)
                - b.mp.getAlpha() / 1.88973 * (a.getCoordinates()[tau] - b.getCoordinates()[tau]) / R * Math.exp(-b.mp.getAlpha() * R / 1.88973);
    }

    public MNDOParams getParams() {

    }

}
