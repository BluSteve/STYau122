package nddoparam.param;

import nddoparam.*;
import nddoparam.mndo.MNDOAtom;
import nddoparam.mndo.MNDOParams;
import scf.AtomHandler;
import scf.Utils;

import java.util.Arrays;

public class ParamGradientRestricted extends ParamGradientAnalytical {

    public ParamGradientRestricted(NDDOSolution s, String kind, int charge, double[] datum, NDDOSolution sExp) {
        super(s, kind, charge, datum, sExp);
        e = new ParamErrorFunctionRestricted(s, datum[0]);
        if (datum[1] != 0) e.addDipoleError(datum[1]);
        if (datum[2] != 0) e.addIEError(datum[2]);

        if (this.sExp != null) {
            isExpAvail = true;
            e.createExpGeom(this.sExp.atoms, this.sExp);
            e.addGeomError();
        }
    }

    // gets only HFDerivs
    @Override
    protected void computeADerivs() {
        for (int Z : s.getUniqueZs()) {
            for (int paramNum : s.getNeededParams()[Z]) {
                HFDerivs[Z][paramNum] = ParamDerivative.HfDeriv((NDDOSolutionRestricted) s, Z, paramNum);
                totalDerivs[Z][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[Z][paramNum];
                if (isExpAvail) computeGeomDeriv(Z, paramNum);
            }
        }
    }

    // gets HFDerivs and dipoleDerivs
    @Override
    protected void computeBDerivs() {
        for (int Z : s.getUniqueZs()) {
            computeBatchedDerivs(Z);

            for (int paramNum : s.getNeededParams()[Z]) {
                computeDipoleDeriv(Z, paramNum, true);
                if (isExpAvail) computeGeomDeriv(Z, paramNum);
            }
        }
    }

    // gets HFDerivs and IEDerivs
    @Override
    protected void computeCDerivs() {
        for (int Z : s.getUniqueZs()) {
            computeBatchedDerivs(Z);

            for (int paramNum : s.getNeededParams()[Z]) {
                computeDipoleDeriv(Z, paramNum, false);
                computeIEDeriv(Z, paramNum);
                if (isExpAvail) computeGeomDeriv(Z, paramNum);
            }
        }
    }

    // gets HFDerivs, dipoleDerivs and IEDerivs;
    @Override
    protected void computeDDerivs() {
        for (int Z : s.getUniqueZs()) {
            computeBatchedDerivs(Z);

            for (int paramNum : s.getNeededParams()[Z]) {
                computeDipoleDeriv(Z, paramNum, true);
                computeIEDeriv(Z, paramNum);
                if (isExpAvail) computeGeomDeriv(Z, paramNum);
            }
        }
    }

    @Override
    protected void computeBatchedDerivs(int Z) {
        staticDerivs[Z] = ParamDerivative.MNDOStaticMatrixDeriv((NDDOSolutionRestricted) s, Z);
        xLimited[Z] = ParamDerivative.xarraylimited((NDDOSolutionRestricted) s, staticDerivs[Z][1]);
    }

    @Override
    protected void computeGeomDeriv(int Z, int paramNum) {
        sPrime = new NDDOSolutionRestricted(Utils.perturbAtomParams(sExp.atoms, paramNum, Z), charge); // sExp perturbed put into a solution
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
        totalDerivs[Z][paramNum] += 0.000049 * e.geomGradient * geomDerivs[Z][paramNum];
    }

    // `full` sets whether dipole itself is actually computed
    @Override
    protected void computeDipoleDeriv(int Z, int paramNum, boolean full) {
        // alpha and eisol have to be computed apart from the rest of HFDerivs. I.e. not as a part of dipole.
        if (paramNum == 0 || paramNum == 7) {
            HFDerivs[Z][paramNum] = ParamDerivative.HfDeriv((NDDOSolutionRestricted) s, Z, paramNum);
        }
        if (staticDerivs[Z][0][paramNum] != null || staticDerivs[Z][1][paramNum] != null) {
            HFDerivs[Z][paramNum] = ParamDerivative.MNDOHfDeriv((NDDOSolutionRestricted) s,
                    staticDerivs[Z][0][paramNum], staticDerivs[Z][1][paramNum]);
            densityDerivs[Z][paramNum] = ParamDerivative.densityDerivativeLimited((NDDOSolutionRestricted) s,
                    xLimited[Z][paramNum]);
            if (full) {
                dipoleDerivs[Z][paramNum] = ParamDerivative.MNDODipoleDeriv((NDDOSolutionRestricted) s,
                        densityDerivs[Z][paramNum], Z, paramNum);

                totalDerivs[Z][paramNum] += 800 * (s.dipole - datum[1]) * dipoleDerivs[Z][paramNum];
            }
        }

        totalDerivs[Z][paramNum] += 2 * (s.hf - datum[0]) * HFDerivs[Z][paramNum];
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
            coeffDerivs[Z][paramNum] = ParamDerivative.HOMOcoefficientDerivativeComplementary(xForIE[Z][paramNum],
                    (NDDOSolutionRestricted) s);
            IEDerivs[Z][paramNum] = ParamDerivative.MNDOIEDeriv((NDDOSolutionRestricted) s,
                    coeffDerivs[Z][paramNum], fockDerivs[Z][paramNum]);

            totalDerivs[Z][paramNum] += 200 * (s.homo + datum[2]) * IEDerivs[Z][paramNum];
        }
    }

    private static String csvify(double[] array) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < array.length - 1; i++) {
            sb.append(array[i]).append(',');
        }
        sb.append(array[array.length - 1]);
        return sb.toString();
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
        NDDOSolution expsoln = new NDDOSolutionRestricted(exp, 0);
        double[] datum = new double[]{-17.9, 0, 13.6};
        NDDOGeometryOptimizationRestricted opt = new NDDOGeometryOptimizationRestricted(new NDDOAtom[]{carbon, atom1, atom2, atom3, atom4}, 0);

        ParamGradientRestricted g = new ParamGradientRestricted(opt.s, "a", 0, datum, expsoln);
        g.computeDerivs();
        System.out.println("Test HF (A) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
        System.out.println("Geom Derivs: " + Arrays.deepToString(g.getGeomDerivs()));

        g = new ParamGradientRestricted(opt.s, "b", 0, datum, expsoln);
        g.computeDerivs();
        System.out.println("Test HF (B) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
        System.out.println("Test Dipole (B) Derivs: " + Arrays.deepToString(g.getDipoleDerivs()));

        g = new ParamGradientRestricted(opt.s, "c", 0, datum, expsoln);
        g.computeDerivs();
        System.out.println("Test HF (C) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
        System.out.println("Test IE (C) Derivs: " + Arrays.deepToString(g.getIEDerivs()));

        g = new ParamGradientRestricted(opt.s, "d", 0, datum, expsoln);
        g.computeDerivs();
        System.out.println("Test HF (D) Derivs: " + Arrays.deepToString(g.getHFDerivs()));
        System.out.println("Test Dipole (D) Derivs: " + Arrays.deepToString(g.getDipoleDerivs()));
        System.out.println("Test IE (D) Derivs: " + Arrays.deepToString(g.getIEDerivs()));
        System.out.println(csvify(g.getHFDerivs()[6]));
    }
}
