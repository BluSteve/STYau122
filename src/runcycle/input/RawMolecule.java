package runcycle.input;

import java.util.Arrays;

public class RawMolecule {
    public int index;
    public String name;
    public boolean uhf;
    public int charge, mult, nElectrons;
    public double[] datum;
    public String kind;
    public int[] uniqueZs;
    public RawAtom[] atoms, expGeom;

    @Override
    public String toString() {
        return "RawMolecule{" +
                "index=" + index +
                ", uhf=" + uhf +
                ", charge=" + charge +
                ", mult=" + mult +
                ", datum=" + Arrays.toString(datum) +
                ", atoms=" + Arrays.toString(atoms) +
                ", expGeom=" + Arrays.toString(expGeom) +
                '}';
    }
}
