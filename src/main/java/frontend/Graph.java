package frontend;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import runcycle.RunIterator;
import runcycle.structs.RunInput;
import runcycle.structs.RunOutput;
import tools.Utils;

import java.awt.*;
import java.io.IOException;

public class Graph extends ApplicationFrame {
	public final XYSeriesCollection dataset;
	public final XYSeries errorSeries;

	public Graph(String title) {
		super(title);

		dataset = new XYSeriesCollection();
		errorSeries = new XYSeries("error");
		dataset.addSeries(errorSeries);
		JFreeChart lineChart = ChartFactory.createXYLineChart(title, "Run Number", "Total Error",
						dataset, PlotOrientation.VERTICAL, true, true, false);
		lineChart.getXYPlot().setRangeAxisLocation(AxisLocation.TOP_OR_RIGHT);
		ChartPanel cp = new ChartPanel(lineChart);
		setContentPane(cp);
	}

	public static void main(String[] args) throws InterruptedException {
		Graph graph = new Graph("Total Error against Run Number");
		graph.setPreferredSize(new Dimension(1200, 800));
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
			graph.errorSeries.add(i, ro.ttError);
			i++;
		}
	}
}
