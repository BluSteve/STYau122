package mndoparam.mndo;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class AtomHandler {
    public static AtomProperty[] atoms = new AtomProperty[119];

    public static void populateAtoms() {
        try {
            Scanner s = new Scanner(new File("atom_properties.csv"));
            s.nextLine(); // skips over header
            while (s.hasNext()) {
                String[] strs = s.nextLine().split(",");
                atoms[Integer.parseInt(strs[0])] = new AtomProperty(strs[1], Double.parseDouble(strs[2]), Double.parseDouble(strs[3]));
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
