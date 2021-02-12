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
        return NDDODerivative.gradientUnrestricted(atoms, (NDDOSolutionUnrestricted) s, i, j);
    }

    protected DoubleMatrix[] routine() {
        DoubleMatrix[][] matrices = NDDODerivative.gradientroutine(atoms, (NDDOSolutionUnrestricted) s);

        DoubleMatrix gradient = matrices[0][0];

        DoubleMatrix hessian = NDDOSecondDerivative.hessianroutine(atoms, (NDDOSolutionUnrestricted) s, matrices[1], matrices[2]);

        //DoubleMatrix hessian = DoubleMatrix.eye(gradient.length);

        return new DoubleMatrix[] {gradient, hessian};
    }


}
