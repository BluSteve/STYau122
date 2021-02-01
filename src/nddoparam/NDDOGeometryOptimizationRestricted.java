package nddoparam;

import org.jblas.DoubleMatrix;

public class NDDOGeometryOptimizationRestricted extends NDDOGeometryOptimization {

    @Override
    protected void updateNDDOSolution() {
        s = new NDDOSolutionRestricted(atoms, charge);
    }


    public NDDOGeometryOptimizationRestricted(NDDOAtom[] atoms, int charge) {
        super(atoms, charge, 1);
    }

    protected double derivative(int i, int j) {
        for (int a = 0; a < atoms.length; a++) {
            for (int b = a; b < atoms.length; b++) {
                for (int tau1 = 0; tau1 < 3; tau1++) {
                    for (int tau2 = 0; tau2 < 3; tau2++) {
                        double finite = NDDOSecondDerivative.hessianfinite(atoms, (NDDOSolutionRestricted) s, a, tau1, b, tau2);
                        double analytic = NDDOSecondDerivative.hessian(atoms, (NDDOSolutionRestricted) s, a, tau1, b, tau2);

                        if (Math.abs(finite - analytic) > 1E-4) {
                            System.err.println ("You expected this!");
                            System.err.println (finite + ", " + analytic);
                            System.exit(0);
                        }
                        else {
                            System.err.println (finite + ", " + analytic);
                        }
                    }
                }
            }
        }

        return NDDODerivative.gradient(atoms, (NDDOSolutionRestricted) s, i, j);
    }
}
