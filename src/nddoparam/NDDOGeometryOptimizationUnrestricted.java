package nddoparam;

import org.jblas.DoubleMatrix;

public class NDDOGeometryOptimizationUnrestricted extends NDDOGeometryOptimization {

    @Override
    protected void updateNDDOSolution() {
        s = new NDDOSolutionUnrestricted(atoms, charge, mult);
    }

    public NDDOGeometryOptimizationUnrestricted(NDDOAtom[] atoms, int charge, int mult) {
        super(atoms, charge, mult);
    }

    protected double derivative(int i, int j) {
        return NDDODerivative.gradientUnrestricted(atoms, s, i, j);
    }

    protected DoubleMatrix[] routine() {
        return null;
    }


}
