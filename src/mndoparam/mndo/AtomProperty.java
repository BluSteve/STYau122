package mndoparam.mndo;

public class AtomProperty {
    private static final double HEATCONV = 4.3363E-2;
    private String name;
    private double mass, heat;

    public AtomProperty(String name, double mass, double heat) {
        this.name = name;
        this.mass = mass;
        this.heat = heat;
    }

    public double getHeat() {
        return heat * HEATCONV;
    }

    public void setHeat(double heat) {
        this.heat = heat;
    }

    public double getMass() {
        return mass;
    }

    public void setMass(double mass) {
        this.mass = mass;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
}
