package mndoparam.mndo;

import java.io.Serializable;

public class MNDOParams implements Serializable {
    private double[] params = new double[13];

    public MNDOParams(double alpha, double betas, double betap, double uss, double upp, double zetas, double zetap,
                      double eisol, double gss, double gsp, double hsp, double gpp, double gp2) {
        params[0] = alpha;
        params[1] = betas;
        params[2] = betap;
        params[3] = uss;
        params[4] = upp;
        params[5] = zetas;
        params[6] = zetap;
        params[7] = eisol;
        params[8] = gss;
        params[9] = gsp;
        params[10] = hsp;
        params[11] = gpp;
        params[12] = gp2;
    }

    public MNDOParams(double[] params) {
        this.params = params.clone();
    }

    public MNDOParams() {

    }

    public double getAlpha() {
        return params[0];
    }

    public double getBetas() {
        return params[1];
    }

    public double getBetap() {
        return params[2];
    }

    public double getUss() {
        return params[3];
    }

    public double getUpp() {
        return params[4];
    }

    public double getZetas() {
        return params[5];
    }

    public double getZetap() {
        return params[6];
    }

    public double getEisol() {
        return params[7];
    }

    public double getGss() {
        return params[8];
    }

    public double getGsp() {
        return params[9];
    }

    public double getHsp() {
        return params[10];
    }

    public double getGpp() {
        return params[11];
    }

    public double getGp2() {
        return params[12];
    }

    public void modifyParam(int index, double amnt) {
        params[index] += amnt;
    }

    @Override
    public MNDOParams clone() {
        return new MNDOParams(getAlpha(), getBetas(), getBetap(), getUss(), getUpp(), getZetas(), getZetap(), getEisol(), getGss(), getGsp(), getHsp(), getGpp(), getGp2());
    }
}
