package runcycle.input;

import java.util.Arrays;

public class RawInput {
    public String trainingSet;
    public int[] atomTypes;
    public RawMolecule[] molecules;

    @Override
    public String toString() {
        return "RawInput{" +
                "trainingSet='" + trainingSet + '\'' +
                ", atomTypes=" + Arrays.toString(atomTypes) +
                ", molecules=" + Arrays.toString(molecules) +
                '}';
    }
}
