package nddoparam.mndo;

import nddoparam.NDDOParams;

import java.io.Serializable;

public class MNDOParams implements Serializable { // this doesn't extend NDDOParams because although this is a special
                                                  // case where it is exactly the same as NDDOParams, its siblings could be very different.
    public NDDOParams nddoParams;

    public MNDOParams(double alpha, double betas, double betap, double uss, double upp, double zetas, double zetap,
                      double eisol, double gss, double gsp, double hsp, double gpp, double gp2) {
        nddoParams = new NDDOParams(alpha, betas, betap, uss, upp, zetas, zetap, eisol, gss, gsp, hsp, gpp, gp2);
    }

    public MNDOParams(double[] params) {
        nddoParams = new NDDOParams(params.clone());
    }

    public MNDOParams() {
    }

    public double getAlpha() {
        return nddoParams.getAlpha();
    }

    public double getBetas() {
        return nddoParams.getBetas();
    }

    public double getBetap() {
        return nddoParams.getBetap();
    }

    public double getUss() {
        return nddoParams.getUss();
    }

    public double getUpp() {
        return nddoParams.getUpp();
    }

    public double getZetas() {
        return nddoParams.getZetas();
    }

    public double getZetap() {
        return nddoParams.getZetap();
    }

    public double getEisol() {
        return nddoParams.getEisol();
    }

    public double getGss() {
        return nddoParams.getGss();
    }

    public double getGsp() {
        return nddoParams.getGsp();
    }

    public double getHsp() {
        return nddoParams.getHsp();
    }

    public double getGpp() {
        return nddoParams.getGpp();
    }

    public double getGp2() {
        return nddoParams.getGp2();
    }

    public void modifyParam(int index, double amnt) {
        nddoParams.modifyParam(index, amnt);
    }

    @Override
    public MNDOParams clone() {
        return new MNDOParams(getAlpha(), getBetas(), getBetap(), getUss(), getUpp(), getZetas(), getZetap(), getEisol(), getGss(), getGsp(), getHsp(), getGpp(), getGp2());
    }
}
