package nddoparam;

import org.apache.commons.lang.time.StopWatch;
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

        DoubleMatrix[] matrices = routine();

        gradient = matrices[0];
        DoubleMatrix B = matrices[1];

        int index;

        double sum = 0;
        for (int i = 0; i < gradient.length; i++) {
            sum += gradient.get(i) * gradient.get(i);
        }
        sum = Math.sqrt(sum);

        DoubleMatrix searchdir = Solve.pinv(B).mmul(gradient).mul(-1 / sum);
        DoubleMatrix oldgrad;

        double energy = refEnergy - 1;
        double scale;
        int count;
        counter = 0;
        int numIt = 0;

        StopWatch sw = new StopWatch();
        sw.start();
        sw.suspend();
        while (mag(gradient) > 0.005) {
            System.out.println("Gradient: " + mag(gradient));

            numIt++;
            scale = 0.1;
            refEnergy = 0;
            double lastEnergyGrad = 0;
            double largestEnergyGrad = 1;
            double energyGrad;
            int numSearch = 0;
            sw.resume();

            while (Math.abs(energy - refEnergy) > 1E-8) {
                refEnergy = energy;
                System.err.println("SCALE" + scale);
                count = 0;

                for (NDDOAtom a : atoms) {
                    for (int i = 0; i < 3; i++) {
                        // get x,y,z of a
                        // positions of atoms are essentially "weights"
                        a.getCoordinates()[i] = Math.round((a.getCoordinates()[i] + scale * searchdir.get(count)) * 1000000000) / 1000000000.0;
                        count++;
                    }
                }

                updateNDDOSolution();

                System.out.println("\nCurrent heat of formation: " + s.hf + "kcal/mol");
                System.out.println("Current HOMO energy: " + s.homo + " eV");
                System.out.println("Current energy: " + s.energy);
                System.out.println("-----------------------------------------------");

                energy = s.energy;

                if (energy > refEnergy) {
                    scale *= -0.5;
                    numSearch = 0;
                }
                else {
                    numSearch++;

//                    if (numSearch >= 5) {
//                        scale *= 2;
//                        numSearch = 0;
//                    }

                    if (lastEnergyGrad == 0) {
                        lastEnergyGrad = energy - refEnergy;
                        largestEnergyGrad = lastEnergyGrad;
                    }
                    else {
                        energyGrad = energy - refEnergy;
                        double x = lastEnergyGrad-energyGrad;
                        double gamma = energyGrad/largestEnergyGrad;
                        double sFactor = (1.1/(1+Math.exp(-3-x))+0.2)*Math.min(1,gamma);
//                        double sFactor = 1.21/(1+Math.exp(5-10*gamma));
//                        double sFactor = (1.2/(1+Math.exp(-3-x)))*(1.1/(1+Math.exp(5-10*gamma)));
                        scale = scale*sFactor;
                        System.out.println(energyGrad + " " + gamma + " " + x + " " + sFactor);
                        lastEnergyGrad = energyGrad;
                        if (energyGrad < largestEnergyGrad)
                            largestEnergyGrad = energyGrad;
                    }
                }
            }
            sw.suspend();
            System.err.println("Time: " + sw.getTime());


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
            DoubleMatrix y = gradient.sub(oldgrad); // difference of gradients

            try {
                B = getb(B, y, searchdir);
                searchdir = Solve.pinv(B).mmul(gradient).mmul(-1);
            } catch (Exception e) {
                System.err.println("Hessian approximation error");
                B = DoubleMatrix.eye(atoms.length * 3);
                searchdir = Solve.pinv(B).mmul(gradient).mmul(-1);
            }

            if (numIt == 5) {
                numIt = 0;
                matrices = routine();
                B = matrices[1];
                searchdir = Solve.pinv(B).mmul(gradient).mul(-1); // TODO is mmul the same as mul when it comes to scalars
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
