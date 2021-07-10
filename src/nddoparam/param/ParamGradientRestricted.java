package nddoparam.param;

import nddoparam.NDDODerivative;
import nddoparam.NDDOGeometryOptimizationRestricted;
import nddoparam.NDDOSolution;
import nddoparam.NDDOSolutionRestricted;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import scf.AtomHandler;
import scf.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class ParamGradientRestricted extends ParamGradientAnalytical {

    public ParamGradientRestricted(NDDOSolution s, String kind, double[] datum, NDDOSolution sExp) {
        super(s, kind, datum, sExp);
        e = new ParamErrorFunctionRestricted(s, datum[0]);
        if (datum[1] != 0) e.addDipoleError(datum[1]);
        if (datum[2] != 0) e.addIEError(datum[2]);

        if (this.sExp != null) {
            isExpAvail = true;
            e.createExpGeom(this.sExp.atoms, this.sExp);
            e.addGeomError();
        }
    }

    @Override
    protected void computeDerivs(String kind) {
        if (kind.equals("b") || kind.equals("c") || kind.equals("d")) computeBatchedDerivs(0, 0);
        for (int Z = 0; Z < s.getUniqueZs().length; Z++) {
            for (int paramNum : s.getNeededParams()[s.getUniqueZs()[Z]]) {
                switch (kind) {
                    case "a":
                        computeHFDeriv(Z, paramNum);
                        break;
                    case "b":
                        computeDipoleDeriv(Z, paramNum, true);
                        break;
                    case "c":
                        computeDipoleDeriv(Z, paramNum, false);
                        computeIEDeriv(Z, paramNum);
                        break;
                    case "d":
                        computeDipoleDeriv(Z, paramNum, true);
                        computeIEDeriv(Z, paramNum);
                        break;
                }
                if (isExpAvail) computeGeomDeriv(Z, paramNum);
            }
        }
    }

    // Compiles all necessary fock matrices into one array before using the Pople algorithm, for faster computation.
    // This function is the only thing that's not computed on a Z, paramNum level.
    // Will not compute anything before the firstZIndex and the firstParamIndex.
    @Override
    protected void computeBatchedDerivs(int firstZIndex, int firstParamIndex) {
        ArrayList<DoubleMatrix> aggregate = new ArrayList<>(s.getUniqueZs().length * NDDOSolution.maxParamNum);
        for (int Z = firstZIndex; Z < s.getUniqueZs().length; Z++) {
            if (Z == firstZIndex)
                staticDerivs[Z] = ParamDerivative.MNDOStaticMatrixDeriv((NDDOSolutionRestricted) s, s.getUniqueZs()[Z],
                        firstParamIndex);
            else staticDerivs[Z] = ParamDerivative.MNDOStaticMatrixDeriv((NDDOSolutionRestricted) s, s.getUniqueZs()[Z],
                    0);
            Collections.addAll(aggregate, staticDerivs[Z][1]);
        }

        DoubleMatrix[] aggregateArray = new DoubleMatrix[aggregate.size()];
        for (int x = 0; x < aggregate.size(); x++) {
            aggregateArray[x] = aggregate.get(x);
        }
        DoubleMatrix[] xLimitedAggregate = ParamDerivative.xArrayLimitedPople((NDDOSolutionRestricted) s,
                aggregateArray);
        int i = 0;
        for (int Z = firstZIndex; Z < s.getUniqueZs().length; Z++) {
            xLimited[Z] = Arrays.copyOfRange(xLimitedAggregate, i * NDDOSolution.maxParamNum,
                    i * NDDOSolution.maxParamNum + NDDOSolution.maxParamNum);
            i++;
        }
    }

    @Override
    protected void computeGeomDeriv(int Z, int paramNum) {
        // TODO sExp.charge or s.charge?
        sPrime = new NDDOSolutionRestricted(Utils.perturbAtomParams(sExp.atoms, paramNum, s.getUniqueZs()[Z]),
                s.charge); // sExp perturbed put into a solution
        double sum = 0;
        double d;
        for (int i = 0; i < sPrime.atoms.length; i++) {
            for (int j = 0; j < 3; j++) {
                d = NDDODerivative.grad((NDDOSolutionRestricted) sPrime, i, j);
                sum += d * d;
            }
        }
        double geomGradient = 627.5 * Math.sqrt(sum);
        geomDerivs[Z][paramNum] = 1 / LAMBDA * (geomGradient - e.geomGradient);
        gradients[Z][paramNum] += 0.000049 * e.geomGradient * geomDerivs[Z][paramNum];
    }

    @Override
    protected void computeHFDeriv(int Z, int paramNum) {
        HFDerivs[Z][paramNum] = ParamDerivative.HFDeriv((NDDOSolutionRestricted) s, s.getUniqueZs()[Z], paramNum);
        gradients[Z][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[Z][paramNum];
    }

    // `full` sets whether dipole itself is actually computed
    @Override
    protected void computeDipoleDeriv(int Z, int paramNum, boolean full) {
        // alpha and eisol have to be computed apart from the rest of HFDerivs. I.e. not as a part of dipole.
        if (paramNum == 0 || paramNum == 7) {
            HFDerivs[Z][paramNum] = ParamDerivative.HFDeriv((NDDOSolutionRestricted) s, s.getUniqueZs()[Z], paramNum);
        }
        if (staticDerivs[Z][0][paramNum] != null || staticDerivs[Z][1][paramNum] != null) {
            HFDerivs[Z][paramNum] = ParamDerivative.MNDOHFDeriv((NDDOSolutionRestricted) s,
                    staticDerivs[Z][0][paramNum], staticDerivs[Z][1][paramNum]);
            densityDerivs[Z][paramNum] = ParamDerivative.densityDerivativeLimited((NDDOSolutionRestricted) s,
                    xLimited[Z][paramNum]);
            if (full) {
                dipoleDerivs[Z][paramNum] = ParamDerivative.MNDODipoleDeriv((NDDOSolutionRestricted) s,
                        densityDerivs[Z][paramNum], s.getUniqueZs()[Z], paramNum);

                gradients[Z][paramNum] += 800 * (s.dipole - datum[1]) * dipoleDerivs[Z][paramNum];
            }
        }

        gradients[Z][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[Z][paramNum];
    }

    @Override
    protected void computeIEDeriv(int Z, int paramNum) {
        if (staticDerivs[Z][0][paramNum] != null || staticDerivs[Z][1][paramNum] != null) {
            responseDerivs[Z][paramNum] = ParamDerivative.responseMatrix((NDDOSolutionRestricted)
                    s, densityDerivs[Z][paramNum]);
            fockDerivs[Z][paramNum] = staticDerivs[Z][1][paramNum].add(responseDerivs[Z][paramNum]);
            xComplementary[Z][paramNum] = ParamDerivative.xArrayComplementary((NDDOSolutionRestricted) s,
                    fockDerivs[Z][paramNum]);
            xForIE[Z][paramNum] = ParamDerivative.xarrayForIE((NDDOSolutionRestricted) s,
                    xLimited[Z][paramNum], xComplementary[Z][paramNum]);
            coeffDerivs[Z][paramNum] = ParamDerivative.HOMOCoefficientDerivativeComplementary(xForIE[Z][paramNum],
                    (NDDOSolutionRestricted) s);
            IEDerivs[Z][paramNum] = -ParamDerivative.MNDOIEDeriv((NDDOSolutionRestricted) s,
                    coeffDerivs[Z][paramNum], fockDerivs[Z][paramNum]);

            gradients[Z][paramNum] += 200 * (s.homo + datum[2]) * IEDerivs[Z][paramNum];
        }
    }

    public static void main(String[] args) {
        AtomHandler.populateAtoms();
        MNDOParams h = new MNDOParams(2.92397599125172,
                -6.222578482830868, 0.0, -12.200235077462583, 0.0,
                1.0693232546199132, 0.0, -13.00142320543855, 12.848,
                0.0, 0.0, 0.0, 0.0);
        MNDOParams c = new MNDOParams(2.5572499654157435, -18.854021376560777, -8.377666892780198, -52.57072065877964, -39.05266019981942, 1.838438013363027, 1.805140784089995, -120.60738371097112, 12.23, 11.47, 2.43, 11.08, 9.84);
        MNDOAtom atom1 = new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{0.635 * Utils.bohr, 0.639 * Utils.bohr, 0.635 * Utils.bohr}, h);
        MNDOAtom atom2 = new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{0.635 * Utils.bohr, -0.635 * Utils.bohr, -0.639 * Utils.bohr}, h);
        MNDOAtom atom3 = new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{-0.639 * Utils.bohr, -0.635 * Utils.bohr, 0.635 * Utils.bohr}, h);
        MNDOAtom atom4 = new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{-0.639 * Utils.bohr, 0.639 * Utils.bohr, -0.639 * Utils.bohr}, h);
        MNDOAtom carbon = new MNDOAtom(AtomHandler.atomsMap.get("C"), new double[]{-0.0021 * Utils.bohr, 0.0021 * Utils.bohr, -0.0021 * Utils.bohr}, c);

        MNDOAtom[] exp = new MNDOAtom[]{
                new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{0.6304 * Utils.bohr, 0.6304 * Utils.bohr, 0.6304 * Utils.bohr}, h),
                new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{0.6304 * Utils.bohr, -0.6304 * Utils.bohr, -0.6304 * Utils.bohr}, h),
                new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{-0.6304 * Utils.bohr, -0.6304 * Utils.bohr, 0.6304 * Utils.bohr}, h),
                new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{-0.6304 * Utils.bohr, 0.6304 * Utils.bohr, -0.6304 * Utils.bohr}, h),
                new MNDOAtom(AtomHandler.atomsMap.get("C"), new double[]{0, 0, 0}, c)};
        MNDOAtom[] exp1 = new MNDOAtom[]{
                new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{0.6304 * Utils.bohr, 0.6304 * Utils.bohr, 0.6304 * Utils.bohr}, h),
                new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{0.6304 * Utils.bohr, -0.6304 * Utils.bohr, -0.6304 * Utils.bohr}, h),
                new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{-0.6304 * Utils.bohr, -0.6304 * Utils.bohr, 0.6304 * Utils.bohr}, h),
                new MNDOAtom(AtomHandler.atomsMap.get("H"), new double[]{-0.6304 * Utils.bohr, 0.6304 * Utils.bohr, -0.6304 * Utils.bohr}, h),
                new MNDOAtom(AtomHandler.atomsMap.get("C"), new double[]{0, 0, 0}, c)};
        NDDOSolution expsoln = new NDDOSolutionRestricted(exp, 0);
        double[] datum = new double[]{-17.9, 0, 13.6};

        NDDOGeometryOptimizationRestricted opt = new NDDOGeometryOptimizationRestricted(exp1, 0);

        ParamGradientRestricted g = new ParamGradientRestricted(opt.s, "a", datum, expsoln);
        g.computeDerivs();
        System.out.println("Test HF (A) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
        System.out.println("Geom Derivs: " + Arrays.deepToString(g.getGeomDerivs()));

        g = new ParamGradientRestricted(opt.s, "b", datum, expsoln);
        g.computeDerivs();
        System.out.println("Test HF (B) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
        System.out.println("Test Dipole (B) Derivs: " + Arrays.deepToString(g.getDipoleDerivs()));

        g = new ParamGradientRestricted(opt.s, "c", datum, expsoln);
        g.computeDerivs();
        System.out.println("Test HF (C) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
        System.out.println("Test IE (C) Derivs: " + Arrays.deepToString(g.getIEDerivs()));

        StopWatch sw = new StopWatch();
        sw.start();
        g = new ParamGradientRestricted(opt.s, "d", datum, expsoln);
        g.computeDerivs();
        sw.stop();
        System.out.println("D Time taken: " + sw.getTime());
        System.out.println("Test HF (D) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
        System.out.println("Test Dipole (D) Derivs: " + Arrays.deepToString(g.getDipoleDerivs()));
        System.out.println("Test IE (D) Derivs: " + Arrays.deepToString(g.getIEDerivs()));
        System.out.println("Test Total (D) Derivs: " + Arrays.deepToString(g.getGradients()));

        sw.reset();
        sw.start();
        ParamHessianRestricted hessian = new ParamHessianRestricted(opt.s, "b", datum, expsoln);
        hessian.computeHessian();
        sw.stop();
        System.out.println("Hessian time taken: " + sw.getTime());
        System.out.println("Test Hessian: " + Arrays.deepToString(hessian.getHessian()));
    }
}
