package nddoparam;

import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.jblas.Solve;

public abstract class GeometryOptimization {
	public double IE, dipole, heat;
	public int charge;
	public Solution s;
	protected NDDOAtom[] atoms;
	protected double refEnergy;
	protected int counter, mult;
	protected DoubleMatrix gradient;

	public GeometryOptimization(NDDOAtom[] atoms, int charge, int mult) {
		this.atoms = atoms;
		this.charge = charge;
		this.mult = mult;
		updateNDDOSolution();
		refEnergy = s.energy;
		System.out.println(
				"\n" + s.getMoleculeName() + " Current heat of formation: " + s.hf +
						"kcal/mol");
		System.out
				.println(s.getMoleculeName() + " Current HOMO energy: " + s.homo + " " +
						"eV");
		System.out.println("-----------------------------------------------");

		DoubleMatrix[] matrices = routine();

		gradient = matrices[0];
		DoubleMatrix B = matrices[1];

		int index;


		DoubleMatrix[] ms = Eigen.symmetricEigenvectors(B);

		DoubleMatrix h = ms[1].diag();

		DoubleMatrix U = ms[0];

		int counter = 0;

		for (int i = 0; i < h.length; i++) {
			if (Math.abs(h.get(i)) > 1E-5) {
				counter++;
			}
		}


		double lambda = lambda(h, U.transpose().mmul(gradient), h.rows - counter);

//        System.err.println("lambda: " + lambda);


		DoubleMatrix searchdir =
				Solve.pinv(B.sub(DoubleMatrix.eye(B.rows).mmul(lambda))).mmul(gradient)
						.mul(-1);

		if (mag(searchdir) > 0.3) {
			searchdir = searchdir.mmul(0.3 / mag(searchdir));
		}

		DoubleMatrix oldgrad;

		double energy = refEnergy - 1;
		int count;
		int numIt = 0;

		StopWatch sw = new StopWatch();
		sw.start();
		sw.suspend();
		while (mag(gradient) > 0.001) {
			System.out.println("Gradient: " + mag(gradient));

			numIt++;
			refEnergy = 0;
			sw.resume();

			refEnergy = energy;
			count = 0;

			for (NDDOAtom a : atoms) {
				for (int i = 0; i < 3; i++) {
					// get x,y,z of a
					// positions of atoms are essentially "weights"
					a.getCoordinates()[i] = Math.round(
							(a.getCoordinates()[i] + searchdir.get(count)) * 1000000000) /
							1000000000.0;
					count++;
				}
			}

			updateNDDOSolution();

			System.out.println("\nCurrent heat of formation: " + s.hf + "kcal/mol");
			System.out.println("Current HOMO energy: " + s.homo + " eV");
			System.out.println("Current energy: " + s.energy);
			System.out.println("-----------------------------------------------");

			energy = s.energy;
			sw.suspend();
//            System.err.println("Time: " + sw.getTime());


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

			DoubleMatrix y = gradient.sub(oldgrad); // difference of gradients

			try {
				B = getb(B, y, searchdir);
			} catch (Exception e) {
				System.err.println("Hessian approximation error");
				B = DoubleMatrix.eye(atoms.length * 3);
			}

			if (numIt == 20) {
				numIt = 0;
				matrices = routine();
				B = matrices[1]; // TODO is mmul the same as mul when it comes to scalars
			}

			ms = Eigen.symmetricEigenvectors(B);

			h = ms[1].diag();

			U = ms[0];


			lambda = lambda(h, U.transpose().mmul(gradient), h.rows - counter);

//            System.err.println("lambda: " + lambda);

			searchdir = Solve.pinv(B.sub(DoubleMatrix.eye(B.rows).mmul(lambda)))
					.mmul(gradient).mul(-1);

			if (mag(searchdir) > 0.3) {
				searchdir = searchdir.mmul(0.3 / mag(searchdir));
			}
		}
		System.out.println("FINAL:");

		updateNDDOSolution();

		System.out.println("\nHeat of formation: " + s.hf + "kcal/mol");
		System.out.println("HOMO energy: " + s.homo + " eV");

		this.dipole = s.dipole;
		this.heat = s.hf;
		this.IE = -s.homo;
	}

	private static double lambda(DoubleMatrix h, DoubleMatrix g, int count) {

		if (count == h.length) {
			return 0;
		}

		double initialguess = h.get(count) - 3;

		double newguess = initialguess + 2;

		while (Math.abs(initialguess - newguess) > 1E-7) {

			initialguess = newguess;

			double f = -initialguess;

			double fprime = -1;

			for (int i = 0; i < h.rows; i++) {
				f += g.get(i) * g.get(i) / (initialguess - h.get(i));
				fprime -= g.get(i) * g.get(i) /
						((initialguess - h.get(i)) * (initialguess - h.get(i)));
			}

			newguess = initialguess - f / fprime;
		}

		return newguess;
	}

	protected abstract void updateNDDOSolution();

	private DoubleMatrix getb(DoubleMatrix B, DoubleMatrix y, DoubleMatrix searchdir) {
		double a = 1 / y.transpose().mmul(searchdir).get(0);
		double b = searchdir.transpose().mmul(B).mmul(searchdir).get(0);
		DoubleMatrix m2 =
				B.mmul(searchdir).mmul(searchdir.transpose()).mmul(B.transpose()).mmul(b);
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