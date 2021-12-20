package nddo;

public interface NDDOOrbitalMethods<T extends NDDOOrbital> {
	double OneCenterERI(T a, T b, T c, T d);

	double G(T a, T b, T c, T d);

	double Ggd(T a, T b, T c, T d, int tau);

	double Gg2d(T a, T b, T c, T d, int tau1, int tau2);

	double Gpd(T a, T b, T c, T d, int num, int type);

	double Gp2d(T a, T b, T c, T d, int num1, int type1, int num2, int type2);

	double Gpgd(T a, T b, T c, T d, int num, int type, int tau);

	double Gp2gd(T a, T b, T c, T d, int num1, int type1, int num2, int type2, int tau);

	double H(T a, T b);

	double Hgd(T a, T b, int tau);

	double Hg2d(T a, T b, int tau1, int tau2);

	double Hzetapd(T a, T b, int num, int type);

	double Hbetapd(T a, T b, int num);

	double Hzetazetap2d(T a, T b, int num1, int type1, int num2, int type2);

	double Hbetazetap2d(T a, T b, int num1, int num2, int type2);

	double Hzetapgd(T a, T b, int num, int type, int tau);

	double Hbetapgd(T a, T b, int num, int tau);

	double Hzetazetap2gd(T a, T b, int num1, int type1, int num2, int type2, int tau);

	double Hbetazetap2gd(T a, T b, int num1, int num2, int type2, int tau);
}
