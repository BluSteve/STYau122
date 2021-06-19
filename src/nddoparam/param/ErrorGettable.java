package nddoparam.param;

public interface ErrorGettable {
    ParamErrorFunction getE();
    void constructErrors(double refHeat);
}
