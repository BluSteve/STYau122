package nddo;

public interface NDDOOrbital {
	NDDOAtom getAtom();

	double U();

	double beta();

	double[] decomposition(double[] point1, double[] point2);

	double[] decomposition2(double[] point1, double[] point2);

	NDDO6G[] orbitalArray();
}
