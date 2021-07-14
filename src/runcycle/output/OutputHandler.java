package runcycle.output;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import runcycle.MoleculeRun;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class OutputHandler {
    public static MoleculeOutput toMoleculeOutput(MoleculeRun result) {
        MoleculeOutput mo = new MoleculeOutput();
        mo.index = result.getRawMolecule().index;
        mo.name = result.getRawMolecule().name;
        mo.time = result.getTime();
        mo.datum = result.getDatum();
        mo.hessian = result.getH().getHessian();

        mo.hf = result.getG().getS().hf;
        mo.dipole = result.getG().getS().dipole;
        mo.ie = -result.getG().getS().homo;
        mo.geomGradient = result.getG().getE().geomGradient;
        mo.total = result.getG().getE().getTotalError();

        ParamGradientOutput pgo = new ParamGradientOutput();
        pgo.hf = result.getG().getHFDerivs();
        pgo.dipole = result.getG().getDipoleDerivs();
        pgo.ie = result.getG().getIEDerivs();
        pgo.geom = result.getG().getGeomDerivs();
        pgo.total = result.getG().getTotalGradients();

        mo.gradient = pgo;
        return mo;
    }

    public static void output(MoleculeOutput[] mos) {
        GsonBuilder builder = new GsonBuilder();
        builder.setPrettyPrinting();
        Gson gson = builder.create();
        try {
            FileWriter fw = new FileWriter("output.json");
            gson.toJson(mos, fw);
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
