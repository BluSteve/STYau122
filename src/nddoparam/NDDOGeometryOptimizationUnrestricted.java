package nddoparam;

import org.jblas.DoubleMatrix;

public class NDDOGeometryOptimizationUnrestricted extends NDDOGeometryOptimization {
    private DoubleMatrix alphaDensity, betaDensity;

    @Override
    protected void updateNDDOSolution() {
        s = new NDDOSolutionUnrestricted(atoms, charge, mult);
    }

    @Override
    protected void updateMatrices() {
        alphaDensity = s.alphaDensity();
        betaDensity = s.betaDensity();
    }

    public NDDOGeometryOptimizationUnrestricted(NDDOAtom[] atoms, int charge, int mult) {
        super(atoms, charge, mult);
    }

    protected double derivative(int i, int j) {
        return NDDODerivative.gradientUnrestricted(atoms, alphaDensity, betaDensity, i, j);
    }


}
