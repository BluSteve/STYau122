package nddo;

public interface NDDOOrbital {
	NDDOAtom getAtom();

	double U();

	int getL();

	NDDOOrbital[] orbitalArray();
}
