package mndoparam.mndo;

import org.jblas.DoubleMatrix;

public class MNDOGeometryOptimizationRestricted extends MNDOGeometryOptimization {
    private DoubleMatrix densityMatrix;

    @Override
    protected void updateMNDOSolution() {
        s = new MNDOSolutionRestricted(atoms, charge);
    }

    @Override
    protected void updateMatrices() {
        densityMatrix = s.densityMatrix();
    }

    public MNDOGeometryOptimizationRestricted(MNDOAtom[] atoms, int charge) {
        super(atoms, charge, 1);
    }

    protected double derivative(int i, int j) {
        return MNDODerivative.gradient(atoms, densityMatrix, i, j);
    }

}
