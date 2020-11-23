package optimize;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import org.jblas.*;

import datum.*;

public class ParamOptimizer {
	
	private ArrayList <ReferenceData> datum;
	
	private double value;
	
	public double[] changes;

	public ParamOptimizer() {
		
		this.datum = new ArrayList <ReferenceData>();
		
	}
	
	public void addData (ReferenceData data) {
		this.datum.add(data);
		
		this.value += data.getValue();
	}
	
	public void optimize (DoubleMatrix B, DoubleMatrix gradient) throws Exception {
		
		PrintWriter pw = new PrintWriter(new FileOutputStream(new File ("mndooutput.txt"), true));
		
		DoubleMatrix searchdir = Solve.pinv(B).mmul(gradient);
		
		//DoubleMatrix searchdir = gradient;
		
		double sum = 0;
		
		for (int i = 0; i < searchdir.rows; i++) {
			sum += searchdir.get(i) * searchdir.get(i);
		}
		
		sum = Math.sqrt(sum);
		
		searchdir = searchdir.mmul(1/sum);
		
		double k = -0.001;
		
		double lambda = 0;
		
		double val = 0;
		
		this.changes = new double[searchdir.rows];
		
		int count = 0;
		
		while (Math.abs(val - value) > 1E-6 && Math.abs(lambda) <= 0.05) {
			
			count++;
			
			lambda += k;
			
			val = value;
			
			changes = searchdir.dup().mmul(lambda).toArray();
			
			value = 0;
			
			
			for (ReferenceData d: datum) {
				d.update(changes);
				value += d.getValue();
			}
			
			if (value > val) {
				k *= -0.5;
			}
			
			
			pw.println ("evaluating: " + lambda + ", " +  value + ", " + Arrays.toString(changes));
		}
		
		pw.println ("FINAL: " + value + ", " + Arrays.toString(changes));
		
		pw.close();
		
	}

}
