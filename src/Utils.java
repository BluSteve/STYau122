public class Utils {
    public static double[] toDoubles(String[] strs) {
        double[] doubles = new double[strs.length];
        for (int i = 0; i < strs.length; i++) {
            doubles[i] = Double.parseDouble(strs[i]);
        }
        return doubles;
    }
}
