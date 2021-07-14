package runcycle.input;

import java.util.Arrays;

public class RawInput {
    public String trainingSet;
    public RawMolecule[] molecules;

    @Override
    public String toString() {
        return "RawInput{" +
                "trainingSet='" + trainingSet + '\'' +
                ", molecules=" + Arrays.toString(molecules) +
                '}';
    }
}
