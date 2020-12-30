package nddoparam;

import org.jblas.DoubleMatrix;

public class NDDOGeometryOptimizationRestricted extends NDDOGeometryOptimization {
    private DoubleMatrix densityMatrix;

    @Override
    protected void updateNDDOSolution() {
        s = new NDDOSolutionRestricted(atoms, charge);
    }

    @Override
    protected void updateMatrices() {
        densityMatrix = s.densityMatrix();
    }

    public NDDOGeometryOptimizationRestricted(NDDOAtom[] atoms, int charge) {
        super(atoms, charge, 1);
    }

    protected double derivative(int i, int j) {
        return NDDODerivative.gradient(atoms, densityMatrix, i, j);
    }
}
