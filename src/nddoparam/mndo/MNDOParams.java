package nddoparam.mndo;

import nddoparam.NDDOParams;

import java.io.Serializable;
import java.util.Arrays;

public class MNDOParams extends NDDOParams{

    public MNDOParams(double alpha, double betas, double betap, double uss, double upp, double zetas, double zetap,
                      double eisol, double gss, double gsp, double hsp, double gpp, double gp2) {
        super(alpha, betas, betap, uss, upp, zetas, zetap, eisol, gss, gsp, hsp, gpp, gp2);
    }
    public MNDOParams(double[] params) {
        super(params.clone());
    }

    public MNDOParams() {}

    @Override
    public MNDOParams clone() {
        return new MNDOParams(getAlpha(), getBetas(), getBetap(), getUss(), getUpp(), getZetas(), getZetap(), getEisol(), getGss(), getGsp(), getHsp(), getGpp(), getGp2());
    }

    @Override
    public String toString() {
        return "MNDOParams{" +
                "params=" + Arrays.toString(params) +
                '}';
    }
}
