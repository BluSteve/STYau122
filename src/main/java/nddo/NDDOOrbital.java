package nddo;

public interface NDDOOrbital<T extends NDDOAtom, S extends NDDOOrbital> {
	T getAtom();

	double U();

	int getL();

	S[] orbitalArray();
}
