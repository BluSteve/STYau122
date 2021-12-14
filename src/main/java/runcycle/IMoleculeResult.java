package runcycle;

import runcycle.structs.RunnableMolecule;

public interface IMoleculeResult {
	RunnableMolecule getUpdatedRm();
	long getTime();
	boolean isExpAvail();

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

	double[][] getHessian();
}
