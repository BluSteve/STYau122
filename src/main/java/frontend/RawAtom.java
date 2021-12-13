package frontend;

import java.util.Arrays;

/* RawAtom doesn't need AtomProperties because AtomProperties
	is constant and the user doesn't need to supply it.
	This is for input.
*/
public class RawAtom {
	public int Z;
	public double[] coords;

	@Override
	public String toString() {
		return "RawAtom{" +
				", Z=" + Z +
				", coords=" + Arrays.toString(coords) +
				'}';
	}
}
