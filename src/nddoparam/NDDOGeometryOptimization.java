package nddoparam;

import org.jblas.DoubleMatrix;
import org.jblas.Solve;
import scf.Utils;

public abstract class NDDOGeometryOptimization {
    protected NDDOAtom[] atoms;
    protected double refEnergy;
    public double IE, dipole, heat;
    public int charge;
    protected int counter, mult;
    public NDDOSolution s;
    protected DoubleMatrix gradient;

    public NDDOGeometryOptimization(NDDOAtom[] atoms, int charge, int mult) {
        this.atoms = atoms;
        this.charge = charge;
        this.mult = mult;
        updateNDDOSolution();
        refEnergy = s.energy;
        System.out.println("\n" + s.getMoleculeName() + " Current heat of formation: " + s.hf + "kcal/mol");
        System.out.println(s.getMoleculeName() + " Current HOMO energy: " + s.homo + " eV");
        System.out.println("-----------------------------------------------");

        DoubleMatrix B = DoubleMatrix.eye(atoms.length * 3);

        if (mult == 1) {

            DoubleMatrix[] matrices = routine();

            gradient = matrices[0];
            B = matrices[1];

        }
        else {
            gradient = new DoubleMatrix(atoms.length * 3, 1);

            int index = 0;
            for (int a = 0; a < atoms.length; a++) {
                for (int i = 0; i < 3; i++) {
                    gradient.put(index, 0, derivative(a, i));
                    index++;
                }
            }




        }

        int index;


        double sum = 0;
        for (int i = 0; i < gradient.length; i++) {
            sum += gradient.get(i) * gradient.get(i);
        }
        sum = Math.sqrt(sum);

        DoubleMatrix searchdir;
        searchdir = Solve.pinv (B).mmul(gradient).mul(-1 / sum);
        DoubleMatrix oldgrad;

        double energy = refEnergy - 1;
        double scale;
        int count;
        counter = 0;

        while (mag(gradient) > 0.005) {

            System.out.println("Gradient: " + mag(gradient));


            scale = 0.1;


            refEnergy = 0;


            while (Math.abs(energy - refEnergy) > 1E-6) {
                refEnergy = energy;

                count = 0;

                for (NDDOAtom a : atoms) {
                    for (int i = 0; i < 3; i++) {
                        a.getCoordinates()[i] = Math.round((a.getCoordinates()[i] + scale * searchdir.get(count)) * 1000000000) / 1000000000.0;
                        count++;
                    }
                }


                updateNDDOSolution();

                System.out.println("\nCurrent heat of formation: " + s.hf + "kcal/mol");
                System.out.println("Current HOMO energy: " + s.homo + " eV");

                System.out.println("-----------------------------------------------");

                energy = s.energy;


                if (energy > refEnergy) {
                    scale *= -0.5;
                }

            }


            refEnergy = energy;

            oldgrad = gradient.dup();

            gradient = new DoubleMatrix(atoms.length * 3, 1);

            index = 0;

            for (int a = 0; a < atoms.length; a++) {
                for (int i = 0; i < 3; i++) {
                    gradient.put(index, 0, derivative(a, i));
                    index++;
                }
            }


            sum = 0;

            DoubleMatrix y = gradient.sub(oldgrad);

            try {
                B = getb(B, y, searchdir);


                searchdir = Solve.pinv(B).mmul(gradient).mmul(-1);
            } catch (Exception e) {
                B = DoubleMatrix.eye(atoms.length * 3);

                searchdir = Solve.pinv(B).mmul(gradient).mmul(-1);

            }

            for (int i = 0; i < gradient.length; i++) {
                sum += searchdir.get(i) * searchdir.get(i);
            }

            sum = Math.sqrt(sum);
            searchdir = searchdir.mmul(0.1 / sum);
            counter++;
        }
        System.out.println("FINAL:");

        updateNDDOSolution();

        System.out.println("\nHeat of formation: " + s.hf + "kcal/mol");
        System.out.println("HOMO energy: " + s.homo + " eV");

        this.dipole = s.dipole;
        this.heat = s.hf;
        this.IE = -s.homo;
    }


    protected abstract void updateNDDOSolution();

    private DoubleMatrix getb(DoubleMatrix B, DoubleMatrix y, DoubleMatrix searchdir) {

        double a = 1 / y.transpose().mmul(searchdir).get(0);

        double b = searchdir.transpose().mmul(B).mmul(searchdir).get(0);

        DoubleMatrix m2 = B.mmul(searchdir).mmul(searchdir.transpose()).mmul(B.transpose()).mmul(b);

        DoubleMatrix m1 = y.mmul(y.transpose()).mmul(a);


        return B.add(m1).sub(m2);
    }

    private double mag(DoubleMatrix gradient) {

        double sum = 0;
        for (int i = 0; i < gradient.length; i++) {
            sum += gradient.get(i) * gradient.get(i);
        }

        return Math.sqrt(sum);
    }

    public NDDOAtom[] getAtoms() {
        return this.atoms;
    }

    protected abstract double derivative(int i, int j);

    protected abstract DoubleMatrix[] routine();
}
