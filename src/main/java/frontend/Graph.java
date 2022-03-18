package frontend;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.time.DynamicTimeSeriesCollection;
import org.jfree.data.time.Second;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;
import tools.Utils;

import java.io.IOException;

public class Graph extends ApplicationFrame {
	public final DynamicTimeSeriesCollection dataset;

	public Graph(String title) {
		super(title);

		dataset = new DynamicTimeSeriesCollection(1, 2000, new Second());
		dataset.setTimeBase(new Second(0, 0, 0, 23, 1, 2014));
		dataset.addSeries(new float[1], 0, title);
		JFreeChart lineChart = ChartFactory.createTimeSeriesChart(title, "Total Error", "Run Number",
						dataset, true, true, false);
		ChartPanel cp = new ChartPanel(lineChart);
		setContentPane(cp);
	}

	public static void main(String[] args) {
		Graph graph = new Graph("Total Error against Run Number");
		graph.pack();
		RefineryUtilities.centerFrameOnScreen(graph);
		graph.setVisible(true);

		Logger logger = LogManager.getLogger();
		logger.info("Date compiled: {}", Utils.getResource("version.txt"));
		FrontendConfig.init();

		RunInput input;
		try {
			input = TxtIO.readInput("molecules.txt");
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		RunIterator iterator = new RunIterator(input, FrontendConfig.config.num_runs);

		logger.info("Number of runs = {}", FrontendConfig.config.num_runs);

		int i = FrontendConfig.config.starting_run_num;
		while (iterator.hasNext()) {
			RunOutput ro = iterator.next();
			graph.dataset.advanceTime();
			graph.dataset.appendData(new float[]{(float) ro.ttError});
			i++;
		}
	}
}
