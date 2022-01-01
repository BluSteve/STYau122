package runcycle;

import runcycle.structs.RunnableMolecule;

public interface IMoleculeResult {
	RunnableMolecule getUpdatedRm();

	long getTime();

	boolean isExpAvail();

	double getHF(); // todo rename

	double getDipole();

	double getIE();

	double getGeomGradMag();

	double getTotalError();

	double[][] getHFDerivs();

	double[][] getDipoleDerivs();

	double[][] getIEDerivs();

	double[][] getGeomDerivs();

	double[][] getTotalGradients();

	double[][] getHessian();
}
