package nddo;

public interface NDDOOrbitalMethods<T extends NDDOOrbital> {
	double OneCenterERI(T a, T b, T c, T d);
	double getG(T a, T b, T c, T d);

	double getGderiv(T a, T b, T c, T d, int tau);

	double getGderiv2(T a, T b, T c, T d, int tau1, int tau2);

	double getGderiv(T a, T b, T c, T d, int num, int type);

	double beta(T a, T b);

	double betaderiv(T a, T b, int tau);

	double betaderiv2(T a, T b, int tau1, int tau2);

	double betaparamderiv(T a, T b, int num, int type);
	
	double betaparambetaderiv(T a, T b, double sum);
}
