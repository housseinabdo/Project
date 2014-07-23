package classes;
import java.awt.Color;
import java.util.HashMap;
import java.util.Iterator;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;

public class DrowChart extends ApplicationFrame {
	public DrowChart(final String title, HashMap<Double, Double> map,
			String yax, String xAx) {
		super(title);
		final XYDataset dataset = createDataset(map);
		final JFreeChart chart = createChart(dataset, yax, xAx);
		final ChartPanel chartPanel = new ChartPanel(chart);
		chartPanel.setPreferredSize(new java.awt.Dimension(800, 800));
		setContentPane(chartPanel);
	}

	private XYDataset createDataset(HashMap<Double, Double> map) {
		// TODO Auto-generated method stub
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries x = new XYSeries("mass");

		Iterator iterator = map.keySet().iterator();
		while (iterator.hasNext()) {
			Double key = (Double) iterator.next();
			Double value = map.get(key);
			System.out.println(key + " " + value);
			x.add(key, value);
		}
		dataset.addSeries(x);
		return dataset;
	}

	private JFreeChart createChart(final XYDataset dataset, String yAx,
			String xAx) {
		final JFreeChart chart = ChartFactory.createXYLineChart("", yAx, xAx,
				dataset, PlotOrientation.VERTICAL, true, true, true);

		chart.setBackgroundPaint(Color.white);
		final XYPlot plot1 = chart.getXYPlot();
		plot1.setBackgroundPaint(Color.lightGray);
		plot1.setDomainGridlinePaint(Color.white);
		plot1.setRangeGridlinePaint(Color.white);
		System.out.println("finishing");
		return chart;
	}
}
