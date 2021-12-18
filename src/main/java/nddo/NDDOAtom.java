package nddo;

import nddo.structs.AtomProperties;

public interface NDDOAtom<T extends NDDOAtom, S extends NDDOOrbital> { // todo make this interface
	AtomProperties getAtomProperties();

	NDDOParams getParams();

	double[] getCoordinates();

	T withNewParams(NDDOParams np);

	T withNewCoords(double[] coordinates);

	S[] getOrbitals();

	S s();

	double V(S a, S b);

	double Vgd(S a, S b, int tau);

	double Vg2d(S a, S b, int tau1, int tau2);

	double Vpd(S a, S b, int num, int type);

	double p0();

	double p1();

	double p1pd(int type);

	double p1p2d(int type);

	double p2();

	double p2pd(int type);

	double p2p2d(int type);

	double D1();

	double D1pd(int type);

	double D1p2d(int type);

	double D2();

	double D2pd(int type);

	double D2p2d(int type);

	double crf(T b);

	double crfDeriv(T b, int tau);

	double crfDeriv2(T b, int tau1, int tau2);

	double crfAlphaDeriv(T b, int num);

	T copy();
}
