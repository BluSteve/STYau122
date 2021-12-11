package runcycle.input;

import java.util.Arrays;

/* RawAtom doesn't need AtomProperties because AtomProperties
	is constant and the user doesn't need to supply it.
	This is for input.
*/
public class RawAtom {
	public String name;
	public int Z;
	public double[] coords = new double[3];

	@Override
	public String toString() {
		return "RawAtom{" +
				"name='" + name + '\'' +
				", Z=" + Z +
				", coords=" + Arrays.toString(coords) +
				'}';
	}
}
