package mndoparam.param;

import mndoparam.mndo.MNDOAtom;
import mndoparam.mndo.MNDODerivative;
import mndoparam.mndo.MNDOSolution;
import scf.GTO;

import java.util.ArrayList;

public class MNDOParamErrorFunction {

    private double HeatError;

    private double DipoleError;

    private double IEError;

    private double GeomError;

    public double gradient;

    private ArrayList<Double> bonderrors;

    private ArrayList<Double> angleerrors;

    public ArrayList<Double> bondderivatives, anglederivatives;

    public ArrayList<Double> bonds, angles;

    public MNDOSolution soln;

    private MNDOAtom[] atoms;

    private MNDOSolution expsoln;

    public MNDOAtom[] expatoms;


    public MNDOParamErrorFunction(MNDOAtom[] atoms, MNDOSolution soln, double refHeat) {

        this.atoms = atoms;

        this.soln = soln;

        this.HeatError = (soln.hf - refHeat) * (soln.hf - refHeat);

        this.DipoleError = 0;

        this.IEError = 0;

        this.angleerrors = new ArrayList<Double>();
        this.bonderrors = new ArrayList<Double>();
        this.bonds = new ArrayList<Double>();
        this.angles = new ArrayList<Double>();
        this.bondderivatives = new ArrayList<Double>();
        this.anglederivatives = new ArrayList<Double>();

        //System.out.println ("Reference heat of formation: " + refHeat);
        //System.out.println ("Calculated heat of formation: " + soln.hf);
    }

    public void AddDipoleError(double refDipole) {
        //System.out.println ("Reference dipole: " + refDipole);
        //System.out.println ("Calculated dipole: " + soln.dipole);
        this.DipoleError = 400 * (soln.dipole - refDipole) * (soln.dipole - refDipole);
    }

    public void AddIEError(double refIE) {
        //System.out.println ("Reference IE: " + refIE);
        //System.out.println ("Calculated IE: " + -soln.homo);
        this.IEError = 100 * (refIE + soln.homo) * (refIE + soln.homo);
    }

    public void createExpGeom(MNDOAtom[] expatoms, MNDOSolution expsoln) {
        this.expatoms = expatoms;
        this.expsoln = expsoln;
    }

    public void AddGeomError() {

        double sum = 0;

        for (int i = 0; i < expatoms.length; i++) {
            for (int j = 0; j < 3; j++) {
                double d = MNDODerivative.gradient(expatoms, expsoln.densitymatrix(), i, j);
                sum += d * d;
            }
        }

        System.err.println("Geometry Obtained");

        this.gradient = 627.5 * Math.sqrt(sum);

        this.GeomError = 0.000049 * 627.5 * 627.5 * sum;

    }

    public void AddBondError(int atom1, int atom2, double ref) {
        double length = GTO.R(atoms[atom1].getCoordinates(), atoms[atom2].getCoordinates()) / 1.88973;
        //System.out.println ("Reference bond length: " + ref);
        //System.out.println ("Calculated bond length: " + length);
        this.bonds.add(length);


        double[] R = new double[]{expatoms[atom2].getCoordinates()[0] - expatoms[atom1].getCoordinates()[0], expatoms[atom2].getCoordinates()[1] - expatoms[atom1].getCoordinates()[1], expatoms[atom2].getCoordinates()[2] - expatoms[atom1].getCoordinates()[2]};

        double dist = GTO.R(expatoms[atom2].getCoordinates(), expatoms[atom1].getCoordinates());

        double deriv = R[0] / dist * MNDODerivative.gradient(expatoms, expsoln.densitymatrix(), atom2, 0)
                + R[1] / dist * MNDODerivative.gradient(expatoms, expsoln.densitymatrix(), atom2, 1)
                + R[2] / dist * MNDODerivative.gradient(expatoms, expsoln.densitymatrix(), atom2, 2);

        deriv = 627.5 * deriv;


        deriv = 1E-13 * Math.round(deriv * 1E13);

        int numrow = MNDODerivative.Hderiv(expatoms, expsoln.densitymatrix(), atom2, 0).rows;

        //System.err.println (MNDOParamDerivative.Hderivative(expatoms, expsoln.densitymatrix(), atom2, 0));

        //System.err.println (MNDOParamDerivative.Gderivative(expatoms, expsoln.densitymatrix(), atom2, 0));


        for (int num = 0; num < numrow; num++) {


            //System.err.println (MNDOParamDerivative.Hderivative(expatoms, expsoln.densitymatrix(), atom2, 1));

            //System.err.println (MNDOParamDerivative.Gderivative(expatoms, expsoln.densitymatrix(), atom2, 1));

            //System.err.println (MNDOParamDerivative.Hderivative(expatoms, expsoln.densitymatrix(), atom2, 2));

            //System.err.println (MNDOParamDerivative.Gderivative(expatoms, expsoln.densitymatrix(), atom2, 2));
        }


        this.bondderivatives.add(deriv);

        this.bonderrors.add(0.49 * (deriv) * (deriv));

    }

    public void AddAngleError(int atom1, int atom2, int atom3, double ref) {
        double[] vector1 = new double[]{expatoms[atom1].getCoordinates()[0] - expatoms[atom2].getCoordinates()[0],
                expatoms[atom1].getCoordinates()[1] - expatoms[atom2].getCoordinates()[1],
                expatoms[atom1].getCoordinates()[2] - expatoms[atom2].getCoordinates()[2]};

        double[] vector2 = new double[]{expatoms[atom3].getCoordinates()[0] - expatoms[atom2].getCoordinates()[0],
                expatoms[atom3].getCoordinates()[1] - expatoms[atom2].getCoordinates()[1],
                expatoms[atom3].getCoordinates()[2] - expatoms[atom2].getCoordinates()[2]};

        double angle = theta(vector1[0], vector1[1], vector1[2], vector2[0], vector2[1], vector2[2]);

        double[] perpendicularvector = normalizedvector(cross(vector1, vector2));


        double[] coeff = normalizedvector(cross(vector2, perpendicularvector));

        double deriv = coeff[0] * MNDODerivative.gradient(expatoms, expsoln.densitymatrix(), atom2, 0)
                + coeff[1] * MNDODerivative.gradient(expatoms, expsoln.densitymatrix(), atom2, 1)
                + coeff[2] * MNDODerivative.gradient(expatoms, expsoln.densitymatrix(), atom2, 2);

        deriv = 627.5 * deriv;

        deriv = 1E-13 * Math.round(deriv * 1E13);

        this.anglederivatives.add(deriv);

        this.angles.add(angle);

        this.angleerrors.add(0.49 * (deriv) * (deriv));
    }

    private static double theta(double x1, double y1, double z1, double x2, double y2, double z2) {
        double sum = x1 * x2 + y1 * y2 + z1 * z2;

        double summ = mag(new double[]{x1, y1, z1}) * mag(new double[]{x2, y2, z2});

        return Math.acos(sum / summ) * 180 / Math.PI;
    }

    private static double[] cross(double[] a, double[] b) {

        return new double[]{a[1] * b[2] - a[1] * b[2], b[0] * a[2] - a[0] * b[2], a[0] * b[1] - b[0] * a[1]};
    }

    private static double[] normalizedvector(double[] v) {
        double val = mag(v);

        return new double[]{v[0] / val, v[1] / val, v[2] / val};
    }


    public double constructErrorFunction() {
        double sum = HeatError + DipoleError + IEError + GeomError;


        for (Double d : angleerrors) {
            sum += d;
        }

        for (Double d : bonderrors) {
            sum += d;
        }


        return sum;
    }

    private static double mag(double[] vector) {

        double sum = 0;
        for (int i = 0; i < vector.length; i++) {
            sum += vector[i] * vector[i];
        }

        return Math.sqrt(sum);
    }

}
