package runcycle.optimize;

import nddo.param.ParamGradient;
import org.ejml.simple.SimpleMatrix;
import runcycle.IMoleculeResult;
import runcycle.optimize.searchdir.ISDFinder;
import runcycle.structs.InputInfo;
import runcycle.structs.LastRunInfo;
import tools.Utils;

import java.util.ArrayList;

public class LLAOptimizerNew extends ParamOptimizer {
	private final ArrayList<ReferenceData> datum;
	private double valueSum;

	public LLAOptimizerNew(IMoleculeResult[] results, InputInfo info, LastRunInfo lri) {
		super(results, info, lri);
		this.datum = new ArrayList<>();
	}


	public double[] optimize(ISDFinder sdFinder) {
		SimpleMatrix searchdir;
		try {
			searchdir = B.pseudoInverse().mult(g);
			SimpleMatrix eigenvalues = Utils.symEigen(B)[1];
			if (logger.isInfoEnabled()) {
				SimpleMatrix eig = eigenvalues.diag().transposei();

				int negCount = 0;
				for (double v : eig.getDDRM().data) if (v < 0) negCount++;

				logger.info("Hessian eigenvalues ({} negative): {}", negCount, eig.toString().strip());
			}
		} catch (IllegalArgumentException e) {
			logger.warn("Hessian pinv failed; using gradient instead.");
			searchdir = g;
		}

		double sum = 0;

		for (int i = 0; i < searchdir.numRows(); i++) {
			sum += searchdir.get(i) * searchdir.get(i);
		}

		sum = Math.sqrt(sum);
		searchdir = searchdir.scale(1 / sum);


		double k = -0.001;
		double lambda = 0;
		double newValueSum = 0;
		double[] changes = new double[searchdir.numRows()];

		while (Math.abs(newValueSum - valueSum) > 1E-8 && Math.abs(lambda) <= 0.5) {
			lambda += k;
			newValueSum = valueSum;
			changes = searchdir.scale(lambda).getDDRM().data;
			valueSum = 0;

			for (ReferenceData d : datum) {
				d.update(changes);
				valueSum += d.getValue();
			}

			if (valueSum > newValueSum) {
				k *= -0.5;
			}
		}

		logger.info("Final lambda: {}", lambda);

		return changes;
	}

	private void addData(ReferenceData data) {
		this.datum.add(data);
		this.valueSum += data.getValue();
	}

	@Override
	protected void processResult(IMoleculeResult result) {
		int[] moleculeATs = result.getUpdatedRm().mats;
		int[][] moleculeNPs = result.getUpdatedRm().mnps;
		boolean isDepad = true;

		double[] datum = result.getUpdatedRm().datum;

		addData(new ReferenceData(datum[0], result.getHf(),
				ParamGradient.combine(result.getHfDerivs(), info.atomTypes, info.neededParams,
						moleculeATs, moleculeNPs, isDepad),
				ReferenceData.HF_WEIGHT));

		if (datum[1] != 0) {
			addData(new ReferenceData(datum[1], result.getDipole(),
					ParamGradient.combine(result.getDipoleDerivs(), info.atomTypes, info.neededParams,
							moleculeATs, moleculeNPs, isDepad),
					ReferenceData.DIPOLE_WEIGHT));
		}

		if (datum[2] != 0) {
			addData(new ReferenceData(datum[2], result.getIE(),
					ParamGradient.combine(result.getIEDerivs(), info.atomTypes, info.neededParams,
							moleculeATs, moleculeNPs, isDepad),
					ReferenceData.IE_WEIGHT));
		}

		if (result.isExpAvail()) {
			addData(new ReferenceData(0, result.getGeomGradMag(),
					ParamGradient.combine(result.getGeomDerivs(), info.atomTypes, info.neededParams,
							moleculeATs, moleculeNPs, isDepad),
					ReferenceData.GEOM_WEIGHT));
		}
	}
}
