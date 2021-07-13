package runcycle.input;

import java.util.Arrays;

public class RawInput {
    String trainingSet;
    RawMolecule[] molecules;

    @Override
    public String toString() {
        return "RawInput{" +
                "trainingSet='" + trainingSet + '\'' +
                ", molecules=" + Arrays.toString(molecules) +
                '}';
    }
}
