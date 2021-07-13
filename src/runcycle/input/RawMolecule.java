package runcycle.input;

import java.util.Arrays;

public class RawMolecule {
    int index;
    boolean uhf;
    int charge, mult;
    RawAtom[] atoms, expGeom;

    @Override
    public String toString() {
        return "RawMolecule{" +
                "index=" + index +
                ", uhf=" + uhf +
                ", charge=" + charge +
                ", mult=" + mult +
                ", atoms=" + Arrays.toString(atoms) +
                ", expGeom=" + Arrays.toString(expGeom) +
                '}';
    }
}
