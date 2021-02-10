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

        return NDDODerivative.gradient(atoms, (NDDOSolutionRestricted) s, i, j);
    }

    protected DoubleMatrix[] routine() {

        DoubleMatrix[][] matrices = NDDODerivative.gradientroutine(atoms, (NDDOSolutionRestricted) s);

        DoubleMatrix gradient = matrices[0][0];

        //DoubleMatrix hessian = NDDOSecondDerivative.hessianroutine(atoms, (NDDOSolutionRestricted) s, matrices[1]);

        DoubleMatrix hessian = DoubleMatrix.eye(gradient.length);

        return new DoubleMatrix[] {gradient, hessian};
    }
}
