package nddoparam.param;

import nddoparam.*;
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

    public ParamGradientRestricted(NDDOSolutionRestricted s, String kind, double[] datum,
                                   NDDOSolutionRestricted sExp, boolean analytical) {
        super(s, kind, datum, sExp, analytical);
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
    protected void constructSPrime(int ZI, int paramNum) {
        sPrime = new NDDOSolutionRestricted(Utils.perturbAtomParams(s.atoms, s.getUniqueZs()[ZI], paramNum), s.charge);
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
    protected void computeHFDeriv(int ZI, int paramNum) {
        if (analytical)
            HFDerivs[ZI][paramNum] = ParamDerivative.HFDeriv((NDDOSolutionRestricted) s, s.getUniqueZs()[ZI], paramNum);
        else
            HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
        totalGradients[ZI][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
    }

    // `full` sets whether dipole itself is actually computed
    // also computes HF
    @Override
    protected void computeDipoleDeriv(int ZI, int paramNum, boolean full) {
        if (analytical) {
            // alpha and eisol have to be computed apart from the rest of HFDerivs. I.e. not as a part of dipole.
            if (paramNum == 0 || paramNum == 7) {
                HFDerivs[ZI][paramNum] = ParamDerivative.HFDeriv((NDDOSolutionRestricted) s, s.getUniqueZs()[ZI], paramNum);
            }
            if (staticDerivs[ZI][0][paramNum] != null || staticDerivs[ZI][1][paramNum] != null) {
                HFDerivs[ZI][paramNum] = ParamDerivative.MNDOHFDeriv((NDDOSolutionRestricted) s,
                        staticDerivs[ZI][0][paramNum], staticDerivs[ZI][1][paramNum]);
                densityDerivs[ZI][paramNum] = ParamDerivative.densityDerivativeLimited((NDDOSolutionRestricted) s,
                        xLimited[ZI][paramNum]);
                if (full) dipoleDerivs[ZI][paramNum] = ParamDerivative.MNDODipoleDeriv((NDDOSolutionRestricted) s,
                        densityDerivs[ZI][paramNum], s.getUniqueZs()[ZI], paramNum);
            }
        } else {
            HFDerivs[ZI][paramNum] = (sPrime.hf - s.hf) / Utils.LAMBDA;
            if (full) dipoleDerivs[ZI][paramNum] = (sPrime.dipole - s.dipole) / Utils.LAMBDA;
        }
        totalGradients[ZI][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[ZI][paramNum];
        if (full) totalGradients[ZI][paramNum] += 800 * (s.dipole - datum[1]) * dipoleDerivs[ZI][paramNum];
    }

    @Override
    protected void computeIEDeriv(int ZI, int paramNum) {
        if (analytical) {
            if (staticDerivs[ZI][0][paramNum] != null || staticDerivs[ZI][1][paramNum] != null) {
                responseDerivs[ZI][paramNum] = ParamDerivative.responseMatrix((NDDOSolutionRestricted)
                        s, densityDerivs[ZI][paramNum]);
                fockDerivs[ZI][paramNum] = staticDerivs[ZI][1][paramNum].add(responseDerivs[ZI][paramNum]);
                xComplementary[ZI][paramNum] = ParamDerivative.xArrayComplementary((NDDOSolutionRestricted) s,
                        fockDerivs[ZI][paramNum]);
                xForIE[ZI][paramNum] = ParamDerivative.xarrayForIE((NDDOSolutionRestricted) s,
                        xLimited[ZI][paramNum], xComplementary[ZI][paramNum]);
                coeffDerivs[ZI][paramNum] = ParamDerivative.HOMOCoefficientDerivativeComplementary(xForIE[ZI][paramNum],
                        (NDDOSolutionRestricted) s);
                IEDerivs[ZI][paramNum] = -ParamDerivative.MNDOIEDeriv((NDDOSolutionRestricted) s,
                        coeffDerivs[ZI][paramNum], fockDerivs[ZI][paramNum]);
            }
        } else {
            IEDerivs[ZI][paramNum] = -(sPrime.homo - s.homo) / Utils.LAMBDA;
        }
        totalGradients[ZI][paramNum] += 200 * -(s.homo + datum[2]) * IEDerivs[ZI][paramNum];
    }

    @Override
    protected void computeGeomDeriv(int Z, int paramNum) {
        // TODO sExp.charge or s.charge?
        sExpPrime = new NDDOSolutionRestricted(Utils.perturbAtomParams(sExp.atoms, s.getUniqueZs()[Z], paramNum), s.charge);
        double sum = 0;
        double d;
        for (int i = 0; i < sExpPrime.atoms.length; i++) {
            for (int j = 0; j < 3; j++) {
                d = NDDODerivative.grad((NDDOSolutionRestricted) sExpPrime, i, j);
                sum += d * d;
            }
        }
        double geomGradient = 627.5 * Math.sqrt(sum);
        geomDerivs[Z][paramNum] = 1 / LAMBDA * (geomGradient - e.geomGradient);
        totalGradients[Z][paramNum] += 0.000098 * e.geomGradient * geomDerivs[Z][paramNum];
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
//        NDDOSolution expsoln = new NDDOSolutionUnrestricted(exp, 0,1 );
        NDDOSolution expsoln = new NDDOSolutionRestricted(exp, 0 );
        double[] datum = new double[]{365.7, 0, 0};
        NDDOGeometryOptimizationRestricted opt = new NDDOGeometryOptimizationRestricted(exp1, 0);

//        NDDOGeometryOptimizationUnrestricted opt = new NDDOGeometryOptimizationUnrestricted(exp1, 0, 1);
        ParamGradientRestricted g;

//        g = new ParamGradientRestricted((NDDOSolutionRestricted) opt.s, "a", datum, (NDDOSolutionRestricted) expsoln, false);
//        g = new ParamGradientUnrestricted2((NDDOSolutionUnrestricted) opt.s, "a", datum, (NDDOSolutionUnrestricted) expsoln, false);
//        g.computeDerivs();
//        System.out.println("Test HF (A) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
//        System.out.println("Geom Derivs: " + Arrays.deepToString(g.getGeomDerivs()));

//        g = new ParamGradientRestricted(opt.s, "b", datum, expsoln, true);
//        g.computeDerivs();
//        System.out.println("Test HF (B) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
//        System.out.println("Test Dipole (B) Derivs: " + Arrays.deepToString(g.getDipoleDerivs()));
//
//        g = new ParamGradientRestricted(opt.s, "c", datum, expsoln, true);
//        g.computeDerivs();
//        System.out.println("Test HF (C) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
//        System.out.println("Test IE (C) Derivs: " + Arrays.deepToString(g.getIEDerivs()));

        StopWatch sw = new StopWatch();
//        double time = 0;
//        for (int x = 0; x < 100; x++) {
//            sw.start();
//            g = new ParamGradientRestricted((NDDOSolutionRestricted)opt.s, "d", datum, (NDDOSolutionRestricted)expsoln, true);
//            g.computeDerivs();
////            ParamHessianRestricted hessian = new ParamHessianRestricted(opt.s, "b", datum, expsoln);
////            hessian.computeHessian();
//            sw.stop();
//            time += sw.getTime();
//            sw.reset();
//        }
        sw.start();
        g = new ParamGradientRestricted((NDDOSolutionRestricted) opt.s, "d", datum, (NDDOSolutionRestricted) expsoln, true);
//        g = new ParamGradientUnrestricted2((NDDOSolutionUnrestricted) opt.s, "d", datum, (NDDOSolutionUnrestricted) expsoln, false);
        g.computeDerivs();

        sw.stop();
        System.out.println("D Time taken: " + sw.getTime());
        System.out.println("Test HF (D) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
        System.out.println("Test Dipole (D) Derivs: " + Arrays.deepToString(g.getDipoleDerivs()));
        System.out.println("Test IE (D) Derivs: " + Arrays.deepToString(g.getIEDerivs()));
        System.out.println("Test Total (D) Derivs: " + Arrays.deepToString(g.getTotalGradients()));

        sw.reset();
        sw.start();
        ParamHessianRestricted hessian = new ParamHessianRestricted((NDDOSolutionRestricted) opt.s, "a", datum, null);
        hessian.computeHessian();
//        ParamHessianUnrestricted2 hessian = new ParamHessianUnrestricted2((NDDOSolutionUnrestricted) opt.s, "b", datum, (NDDOSolutionUnrestricted) expsoln);
//        hessian.computeHessian();
        sw.stop();
        System.out.println("Hessian time taken: " + sw.getTime());
        System.out.println("Test Hessian: " + Arrays.deepToString(hessian.getHessian()));
    }
}
