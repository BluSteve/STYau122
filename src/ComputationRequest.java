import nddoparam.mndo.MNDOAtom;

import java.util.Arrays;

public class ComputationRequest {
    public int index;
    public boolean restricted;
    public MNDOAtom[] atoms;
    public int charge;
    public MNDOAtom[] expgeom;
    public double[] datum;
    public boolean hasHessian;
    public int mult;


    public ComputationRequest(boolean restricted, MNDOAtom[] atoms, int charge, int mult, MNDOAtom[] expgeom, double[] datum,
                              boolean hasHessian, int index) {
        this.restricted = restricted;
        this.atoms = atoms;
        this.charge = charge;
        this.expgeom = expgeom;
        this.datum = datum;
        this.hasHessian = hasHessian;
        this.mult = mult;
        this.index = index;
    }

    @Override
    public String toString() {
        return "ComputationRequest{" + "restricted=" + restricted + ", atoms=" + Arrays
                .toString(atoms) + ", charge=" + charge + ", expgeom=" + Arrays.toString(expgeom) + ", datum=" + Arrays
                .toString(datum) + ", hasHessian=" + hasHessian + ", mult=" + mult + '}';
    }


}