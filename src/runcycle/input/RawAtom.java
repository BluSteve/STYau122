package runcycle.input;

import java.util.Arrays;

public class RawAtom {
    String name;
    double[] coords = new double[3];

    @Override
    public String toString() {
        return "RawAtom{" +
                "name='" + name + '\'' +
                ", coords=" + Arrays.toString(coords) +
                '}';
    }
}
