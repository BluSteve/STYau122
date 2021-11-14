package runcycle;

import runcycle.input.RawMolecule;

public interface MoleculeResult {
	RawMolecule getRm();
	long getTime();
	boolean isExpAvail();
	double[] getDatum();
	double[][] getHessian();

	double getHF();
	double getDipole();
	double getIE();
	double getGeomGradient();
	double getTotalError();

	double[][] getHFDerivs();
	double[][] getDipoleDerivs();
	double[][] getIEDerivs();
	double[][] getGeomDerivs();
	double[][] getTotalGradients();

}
