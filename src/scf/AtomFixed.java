package scf;

public abstract class AtomFixed {//abstract class. Consists of an array of LCGTO objects (ie. the valence orbitals), a coordinate array, an atomic charge Z and a core charge Q.
    protected LCGTO[] orbitals;
    private double[] coordinates;
    protected int shell;
    private int Z, Q;


    public AtomFixed(LCGTO[] orbitals, double[] coordinates, int Z) {
        this.orbitals = orbitals;
        this.coordinates = coordinates;
        this.Z = Z;

        if (Z > 2) {
            this.shell = 2;
            this.Q = Z - 2;
        } else {
            this.shell = 1;
            this.Q = Z;
        }
    }

    public LCGTO[] getOrbitals() {
        return this.orbitals;
    }

    public double[] getCoordinates() {
        return this.coordinates;
    } // only non constant thing

    public int getZ() {
        return this.Z;
    }

    public int getQ() {
        return this.Q;
    }
}
