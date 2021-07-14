package runcycle.input;

import nddoparam.NDDOAtom;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import scf.AtomHandler;

import java.util.Arrays;

public class RawAtom {
    public String name;
    public int Z;
    public double[] coords = new double[3];

    public MNDOAtom toMNDOAtom(MNDOParams mndoParams) {
        return new MNDOAtom(AtomHandler.atoms[Z], coords, mndoParams);
    }

    @Override
    public String toString() {
        return "RawAtom{" +
                "name='" + name + '\'' +
                ", Z=" + Z +
                ", coords=" + Arrays.toString(coords) +
                '}';
    }
}
