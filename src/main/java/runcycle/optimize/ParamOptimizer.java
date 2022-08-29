package runcycle.optimize;

import nddo.param.ParamGradient;
import nddo.param.ParamHessian;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ejml.simple.SimpleMatrix;
import runcycle.IMoleculeResult;
import runcycle.optimize.searchdir.ISDFinder;
import runcycle.structs.InputInfo;
import runcycle.structs.LastRunInfo;

import java.util.Arrays;

public abstract class ParamOptimizer {
	protected static final Logger logger = LogManager.getLogger();
	protected IMoleculeResult[] results;
	protected InputInfo info;
	protected LastRunInfo lri, newLri;
	protected double ttError;
	protected SimpleMatrix B, g;

	public ParamOptimizer(IMoleculeResult[] results, InputInfo info, LastRunInfo lri) {
		this.results = results;
		this.info = info;

		int paramLength = 0; // combined length of all differentiated params
		for (int[] param : info.neededParams) paramLength += param.length;

		ttError = 0;
		double[] ttGradient = new double[paramLength];
		double[][] ttHessian = new double[paramLength][paramLength];

		for (IMoleculeResult result : results) {
			int[] moleculeATs = result.getUpdatedRm().mats;
			int[][] moleculeNPs = result.getUpdatedRm().mnps;
			boolean isDepad = true;

			ttError += result.getTotalError();

			// ttGradient is sum of totalGradients across molecules
			double[] g = ParamGradient.combine(result.getTotalGradients(), info.atomTypes, info.neededParams,
					moleculeATs, moleculeNPs, isDepad);
			for (int i = 0; i < g.length; i++) {
				ttGradient[i] += g[i];
			}

			double[][] h = ParamHessian.padHessian(result.getHessian(), result.getUpdatedRm().mats,
					info.atomTypes, info.neededParams);

			boolean hasNan = false;
			for (int i = 0; i < h.length; i++) {
				for (int j = 0; j < h[0].length; j++) {
					if (Double.isNaN(h[i][j])) {
						hasNan = true;
					}
					else ttHessian[i][j] += h[i][j];
				}
			}

			if (hasNan) {
				logger.warn("NaN in Hessian! {}: \n{}\n{}", result.getUpdatedRm().debugName(), g,
						Arrays.deepToString(h));
			}

			processResult(result);
		}

		logger.info("Total error: {}", ttError);

		SimpleMatrix g = new SimpleMatrix(ttGradient);
		SimpleMatrix B = new SimpleMatrix(ttHessian);
	}

	protected abstract void processResult(IMoleculeResult result);

	public abstract double[] optimize(ISDFinder sdFinder);

	public LastRunInfo getNewLri() {
		return newLri;
	}
}
