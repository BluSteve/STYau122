package nddoparam.mndo;

import org.jblas.DoubleMatrix;

public class MNDOGeometryOptimizationUnrestricted extends MNDOGeometryOptimization {
    private DoubleMatrix alphaDensity, betaDensity;

    @Override
    protected void updateMNDOSolution() {
        s = new MNDOSolutionUnrestricted(atoms, charge, mult);
    }

    @Override
    protected void updateMatrices() {
        alphaDensity = s.alphaDensity();
        betaDensity = s.betaDensity();
    }

    public MNDOGeometryOptimizationUnrestricted(MNDOAtom[] atoms, int charge, int mult) {
        super(atoms, charge, mult);
    }

    protected double derivative(int i, int j) {
        return MNDODerivative.gradientUnrestricted(atoms, alphaDensity, betaDensity, i, j);
    }


}
