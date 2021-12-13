package runcycle;

import runcycle.structs.RunnableMolecule;

public interface MoleculeResult {
	RunnableMolecule getRm();
	long getTime();
	boolean isExpAvail();
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
