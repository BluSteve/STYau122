package scf;

import nddoparam.NDDOAtom;
import nddoparam.NDDOParams;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

public class Utils {
    public static final double lambda = 1E-7;
    public static final double bohr = 1.88973;

    public static double[] toDoubles(String[] strs) {
        double[] doubles = new double[strs.length];
        for (int i = 0; i < strs.length; i++) {
            doubles[i] = Double.parseDouble(strs[i]);
        }
        return doubles;
    }

    public static NDDOAtom[] perturbAtomParams(NDDOAtom[] atoms, int paramNum, int Z) {
        NDDOAtom[] perturbed = new NDDOAtom[atoms.length];
        try {
            Class<? extends NDDOAtom> c = atoms[0].getClass();
            Constructor ctor = c.getDeclaredConstructor(c, atoms[0].getParams().getClass());
            ctor.setAccessible(true);
            Constructor ctor2 = c.getDeclaredConstructor(c);
            ctor2.setAccessible(true);
            for (int i = 0; i < atoms.length; i++) {
                if (atoms[i].getAtomProperties().getZ() == Z) {
                    NDDOParams params = atoms[i].getParams();
                    params.modifyParam(paramNum, Utils.lambda);
                    perturbed[i] = (NDDOAtom) ctor.newInstance(atoms[i], params);
                }
                else {
                    perturbed[i] = (NDDOAtom) ctor2.newInstance(atoms[i]);
                }
            }
        }
        catch (NoSuchMethodException | IllegalAccessException | InstantiationException | InvocationTargetException e) {
            e.printStackTrace();
        }
        return perturbed;
    }

    public static NDDOAtom[] perturbAtomCoords(NDDOAtom[] atoms, int atomnum, int tau) {
        NDDOAtom[] perturbed = new NDDOAtom[atoms.length];
        try {
            Class<? extends NDDOAtom> c = atoms[0].getClass();
            Constructor ctor = c.getDeclaredConstructor(c, atoms[0].getCoordinates().getClass());
            ctor.setAccessible(true);
            Constructor ctor2 = c.getDeclaredConstructor(c);
            ctor2.setAccessible(true);
            for (int i = 0; i < atoms.length; i++) {
                if (i == atomnum) {
                    double[] coords = atoms[i].getCoordinates().clone();
                    coords[tau] = coords[tau] + 1E-7;
                    perturbed[i] = (NDDOAtom) ctor.newInstance(atoms[i], coords);
                }
                else {
                    perturbed[i] = (NDDOAtom) ctor2.newInstance(atoms[i]);
                }
            }
        }
        catch (NoSuchMethodException | IllegalAccessException | InstantiationException | InvocationTargetException e) {
            e.printStackTrace();
        }
        return perturbed;
    }

    public static boolean containsZ(NDDOAtom[] atoms, int Z) {
        boolean result = false;
        for (NDDOAtom atom: atoms) {
            if (atom.getAtomProperties().getZ() == Z) {
                result = true;
                break;
            }
        }
        return result;
    }
    public static int getTrainingSetSize(String trainingSet) {
        int result = 0;
        int removeH = 0;
        if (trainingSet.contains("H")) {
            result += 5;
            removeH = 1;
        }
        result += (trainingSet.length() - removeH) * 8;
        return result;
    }
    /**
     * The confluent hypergeometric function, usually given the symbol 1F1.
     * For information about the function and use of parameters consult the link:
     * http://mathworld.wolfram.com/ConfluentHypergeometricFunctionoftheFirstKind.html
     */
    public static double confluentHypergeometric1F1(double a, double c, double xin) {
        /*
        takes in the numerator parameter of 1F1 (ap), the denominator parameter of 1F1 (cp) and the value of the argument (z)&
    	returns either the value of the rational approximation of 1F1 when it converges to a given tolerance or,if that given
    	tolerance is not met, will ask for different values of ap, cp & z.
    	Also: A and B	will contain the values of the numerator and denominator polynomials, respectively, for all degrees
    	from 0 to n inclusive where n is the maximum degree for which the values of the polynomials are to be calculated
    	*/


        double x = -xin;
        double[] arrayA = new double[1001]; //the numerator of the rational approx.
        double[] arrayB = new double[1001]; //the denominator of the rational approx.
        double[] arrayR = new double[1001]; //the value of the rational approximations.
        double[] arrayD = new double[1001]; //difference in subsequent rational approx.
        int n = 200; //number of iterations,if you change this number must also change condition of the last if statement
        int n1 = n + 1;
        double tolerance = 1E-10;//specified tolerance

        double zero = 0.0;
        double one = 1.0;
        double two = 2.0;
        double three = 3.0;


        // to understand the following code refer to rational approximation of 1F1(ap; cp; -z)

        double ct1 = a * x / c;
        double xn3 = zero;
        double xn1 = two;
        double z2 = x / two;
        double ct2 = z2 / (one + c);
        double xn2 = one;

        arrayA[0] = one;
        arrayB[0] = one;
        arrayB[1] = one + (one + a) * z2 / c;
        arrayA[1] = arrayB[1] - ct1;
        arrayB[2] = one + (two + arrayB[1]) * (two + a) / three * ct2;
        arrayA[2] = arrayB[2] - (one + ct2) * ct1;
        ct1 = three;
        double xn0 = three;
        //for i=3,...,n the values of arrayA[1+i] and arrayB[1+i] are calculated using the recurrence relations below*/

        for (int i = 3; i <= n; i++) {
            //calculation of the multipliers for the recursion
            ct2 = z2 / ct1 / (c + xn1);
            double g1 = one + ct2 * (xn2 - a);
            ct2 = ct2 * (a + xn1) / (c + xn2);
            double g2 = ct2 * ((c - xn1) + (a + xn0) / (ct1 + two) * z2);
            double g3 = ct2 * z2 * z2 / ct1 / (ct1 - two) * (a + xn2) / (c + xn3) * (a - xn2);

            //the recurrance relations for arrayA[i+1] and arrayB[i+1] are as follows

            arrayB[i] = g1 * arrayB[i - 1] + g2 * arrayB[i - 2] + g3 * arrayB[i - 3];
            arrayA[i] = g1 * arrayA[i - 1] + g2 * arrayA[i - 2] + g3 * arrayA[i - 3];

            xn3 = xn2;
            xn2 = xn1;
            xn1 = xn0;
            xn0 = xn0 + one;
            ct1 = ct1 + two;
        }//end for

        arrayD[n1 - 1] = zero;
        for (int j1 = 1; j1 <= n + 1; j1++) {
            arrayR[j1 - 1] = (arrayA[j1 - 1]) / (arrayB[j1 - 1]);//rational approximation of 1f1
            if (j1 > 1) {
                arrayD[j1 - 2] = (arrayR[j1 - 1]) - (arrayR[j1 - 2]);
            }
            if (j1 >= 5 && Math.abs(arrayD[j1 - 2] / arrayR[j1 - 1]) <= tolerance) {
                //checking for convergence of the rational approximation to a given tolerance
                //if that tolerance is met then exit the loop and return the value of the approximation
                return arrayR[j1 - 1];
            }//end if
            //if that tolerance is not met within the given numberof iterations then the program will
            //ask you to check the values entered
            if (j1 == 1000 && Math.abs(arrayD[j1 - 2] / arrayR[j1 - 1]) > tolerance) {
                throw new RuntimeException("please check your the values a, c & xin");
            }
        }//end for

        return arrayR[n];
    }//end main
}
