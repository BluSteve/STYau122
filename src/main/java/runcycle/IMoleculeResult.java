package runcycle;

import runcycle.structs.IEssentialMolecule;

public interface IMoleculeResult {
	IEssentialMolecule getUpdatedRm();

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
