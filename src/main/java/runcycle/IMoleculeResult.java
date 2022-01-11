package runcycle;

import runcycle.structs.RunnableMolecule;

public interface IMoleculeResult {
	RunnableMolecule getUpdatedRm();

	long getTime();

	boolean isExpAvail();

	double getHf();

	double getDipole();

	double getIE();

	double getGeomGradMag();

	double getTotalError();

	double[][] getHfDerivs();

	double[][] getDipoleDerivs();

	double[][] getIEDerivs();

	double[][] getGeomDerivs();

	double[][] getTotalGradients();

	double[][] getHessian();
}
