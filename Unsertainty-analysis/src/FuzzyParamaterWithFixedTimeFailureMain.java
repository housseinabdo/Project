import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.border.EmptyBorder;
import javax.swing.plaf.SliderUI;
import javax.swing.plaf.basic.BasicBorders.RadioButtonBorder;
import javax.swing.text.html.HTMLDocument.Iterator;
import javax.swing.GroupLayout.Alignment;
import javax.swing.*;

import java.awt.event.*;
import java.awt.Color;

import javax.swing.LayoutStyle.ComponentPlacement;

import java.awt.Canvas;
import java.awt.Label;

import javax.swing.JButton;
import javax.swing.JTable;

import fuzzydl.TriangularFuzzyNumber;
import fuzzydl.exception.FuzzyOntologyException;

import javax.swing.SwingConstants;

import java.awt.Rectangle;

import javax.swing.JRadioButton;

import java.awt.Component;

import javax.swing.Box;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;
import org.jfree.ui.RefineryUtilities;

import classes.DrowChart;
import classes.FixedTimeFuzzyODE;
import classes.MonteCarlo;

import com.panayotis.gnuplot.JavaPlot;
import com.panayotis.gnuplot.plot.DataSetPlot;
import com.panayotis.gnuplot.style.PlotStyle;
import com.panayotis.gnuplot.style.Style;

import java.awt.Font;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import classes.*;
public class FuzzyParamaterWithFixedTimeFailureMain extends JFrame {
	static double MAXmassReleased;
	static double MINmassReleased;
	static double MOYmassReleased;
	static double MAXdeg;
	private JPanel contentPane;
	public JTextField renA;
	public JTextField lamda;
	public JTextField renB;
	public JTextField deltah;
	public JLabel lamda_iLabel;
	private JTextField renC;
	private JRadioButton FuzzyRen;
	private JRadioButton NormalRen;
	private ButtonGroup Ren;
	private ButtonGroup mass;
	private JLabel ResultDis;
	private JTextField M;
	public JCheckBox Gaussian;
	public JCheckBox Normal;
	private JTextArea lblMassReleased1;
	private JLabel lblMaxMass;
	/**
	 * @wbp.nonvisual location=-39,279
	 */
	private final JTextField MC = new JTextField();
	private JTextField Mc;
	private JTextField TR;
	private JTextField X;
	private JTextField open;
	private JTextField out;
	private JTextField tbreak;
	private JTextField tmax;
	private JTextField k;
	private JTextField FcValue;
	private JTextField M1;
	private JTextField M2;
	private JTextField Mc1;
	private JTextField Mc2;
	private JTextField TR1;
	private JTextField TR2;
	private JTextField X1;
	private JTextField X2;
	private JTextField open1;
	private JTextField open2;
	private JTextField out1;
	private JTextField out2;
	private JPanel panel_2;
	private JPanel panel_3;
	private JLabel lblLoi;
	private ButtonGroup groupemon;
	private JLabel lblCounter;
	private JLabel lblCounter_1;
	private JScrollPane scr;
	private JTextField txtlSampleNb;
	static HashMap<Double, Double> MassDegre;
	static HashMap<Double, Double> MassDegreFinal;
	static HashMap<Double, Double> DegreeMass; // for represent the degree as a
												// key
	private JTextField textCut;
	private JLabel Cut;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {

		EventQueue.invokeLater(new Runnable() {
			public void run() {

				try {
					FuzzyParamaterWithFixedTimeFailureMain frame = new FuzzyParamaterWithFixedTimeFailureMain();
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	class DistanceThread extends Thread {
		@Override
		public void run() {
			double lamdai, deltaH;

			lamdai = Double.parseDouble(lamda.getText());

			deltaH = Double.parseDouble(deltah.getText());
			MassDegreFinal = new HashMap<Double, Double>();
			if (FuzzyRen.isSelected()) {
				try {
					TriangularFuzzyNumber m = new TriangularFuzzyNumber(
							MINmassReleased, MOYmassReleased, MAXmassReleased);
					ResultDis.setText("");
					TriangularFuzzyNumber rendement = new TriangularFuzzyNumber(
							Double.parseDouble(renA.getText()),
							Double.parseDouble(renB.getText()),
							Double.parseDouble(renC.getText()));
					TriangularFuzzyNumber dist = CalculateDistanceFuzz(
							rendement, lamdai, m, deltaH);

					System.out.println(dist);

					ResultDis.setText("(" + dist.getA() + " , " + dist.getB()
							+ " , " + dist.getC() + ")");
					double ab = 0;
					double hh = Math.pow(10,
							(-(Double.parseDouble(textCut.getText()))));
					for (double h = 0; h < 1; h = h + 0.1) {
						ArrayList<Double> array = new ArrayList<Double>();
						java.util.Iterator<Double> iterator2 = MassDegre
								.keySet().iterator();
						while (iterator2.hasNext()) {
							Double key = (Double) iterator2.next();
							Double value = MassDegre.get(key);

							System.out.println("le valeur de h est" + ab);
							if ((value >= ab) && (value <= ab + hh)) {

								array.add(key);
							}

						}
						ab += hh;

						if (array.size() != 0) {
							Collections.sort(array);
							double min = array.get(0);
							double max = array.get(array.size() - 1);
							double mindeg = MassDegre.get(min);
							double maxdeg = MassDegre.get(max);

							MassDegreFinal.put(min, mindeg);
							MassDegreFinal.put(max, maxdeg);

						}
					}

					final DrowChart demo = new DrowChart("Multi Line Chart",
							MassDegreFinal, "Mass Released", "Degs");
					demo.pack();
					RefineryUtilities.centerFrameOnScreen(demo);
					demo.setVisible(true);

					java.util.Iterator<Double> iterator = MassDegre.keySet()
							.iterator();
					java.util.Iterator<Double> iterator1 = MassDegre.keySet()
							.iterator();

					ArrayList<String> steps = new ArrayList<String>();
					while (iterator.hasNext()) {
						Double key = (Double) iterator.next();
						TriangularFuzzyNumber tri = CalculateDistance(
								rendement, lamdai, key, deltaH);
						steps.add("Mass: " + key + "---> Distance= ("
								+ tri.getA() + " , " + tri.getB() + " , "
								+ tri.getC() + ")");
						Double value = MassDegre.get(key);
						// System.out.println(key + " " + value);
					}
					try {
						PrintWriter writer = new PrintWriter(
								new File(
										"C:/Users/houssein/Desktop/resultsDistance.txt"),
								"UTF-8");
						for (String step : steps) {
							writer.println(step);
						}
						writer.close();

					} catch (Exception e) {
					}
					;

				} catch (FuzzyOntologyException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			if (NormalRen.isSelected()) {
				double m1 = MAXmassReleased;
				double distan = lamdai
						* (Double.parseDouble(renA.getText()) * ((m1 * deltaH) / 4690));
				System.out.println(distan);
				ResultDis.setText("distance = " + distan);

			}
		}

		public TriangularFuzzyNumber CalculateDistanceFuzz(
				TriangularFuzzyNumber rendement, double lamdai,
				TriangularFuzzyNumber m, double deltaH)
				throws FuzzyOntologyException {
			double a = lamdai
					* (rendement.getA() * ((m.getA() * deltaH) / 4690));
			double b = lamdai
					* (rendement.getB() * ((m.getB() * deltaH) / 4690));
			double c = lamdai
					* (rendement.getC() * ((m.getC() * deltaH) / 4690));
			TriangularFuzzyNumber distance = new TriangularFuzzyNumber(a, b, c);
			return distance;
		}

		public TriangularFuzzyNumber CalculateDistance(
				TriangularFuzzyNumber rendement, double lamdai, double m,
				double deltaH) throws FuzzyOntologyException {
			double a = lamdai * (rendement.getA() * ((m * deltaH) / 4690));
			double b = lamdai * (rendement.getB() * ((m * deltaH) / 4690));
			double c = lamdai * (rendement.getC() * ((m * deltaH) / 4690));
			TriangularFuzzyNumber distance = new TriangularFuzzyNumber(a, b, c);
			return distance;
		}
	}

	class MyThread extends Thread {
		double somme = 0, massReleasedValue;

		ArrayList<Double> MassRelacheeList;
		ArrayList<Double> DegreeOfMass;

		public MyThread() {

		}

		@Override
		public void run() {
			lblMaxMass.setText("");
			lblMassReleased1.setText("");
			int Count = Integer.parseInt(txtlSampleNb.getText());
			MassDegre = new HashMap<Double, Double>();
			DegreeMass = new HashMap<Double, Double>();
			MassRelacheeList = new ArrayList<Double>();
			DegreeOfMass = new ArrayList<Double>();
			TriangularFuzzyNumber TrianM = null, TrianMC = null, TrianTR = null, TrianX = null, TrianOpen = null, TrianOut = null;
			try {
				TrianM = new TriangularFuzzyNumber(Double.parseDouble(M
						.getText()), Double.parseDouble(M1.getText()),
						Double.parseDouble(M2.getText()));
				TrianMC = new TriangularFuzzyNumber(Double.parseDouble(Mc
						.getText()), Double.parseDouble(Mc1.getText()),
						Double.parseDouble(Mc2.getText()));
				TrianTR = new TriangularFuzzyNumber(Double.parseDouble(TR
						.getText()), Double.parseDouble(TR1.getText()),
						Double.parseDouble(TR2.getText()));
				TrianX = new TriangularFuzzyNumber(Double.parseDouble(X
						.getText()), Double.parseDouble(X1.getText()),
						Double.parseDouble(X2.getText()));
				TrianOpen = new TriangularFuzzyNumber(Double.parseDouble(open
						.getText()), Double.parseDouble(open1.getText()),
						Double.parseDouble(open2.getText()));
				TrianOut = new TriangularFuzzyNumber(Double.parseDouble(out
						.getText()), Double.parseDouble(out1.getText()),
						Double.parseDouble(out2.getText()));
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (FuzzyOntologyException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			MonteCarlo mon = new MonteCarlo();// for the random variable
			double Tbreakk, TMAX, KK, FFC;
			Tbreakk = Double.parseDouble(tbreak.getText());
			TMAX = Double.parseDouble(tmax.getText());
			KK = Double.parseDouble(k.getText());
			FFC = Double.parseDouble(FcValue.getText());
			int monteCarlo = 1;
			if (Normal.isSelected())
				monteCarlo = 1;
			if (Gaussian.isSelected())
				monteCarlo = 2;

			for (int counter = 0; counter < Count; counter++) {
				lblCounter.setText((counter + 1) + "");
				massReleasedValue = 0;
				somme = 0;
				FirstOrderIntegrator dp853 = new DormandPrince853Integrator(
						6.14e-20, 1.0, 1.0e-10, 1.0e-10);
				FirstOrderDifferentialEquations ode = new FixedTimeFuzzyODE(Tbreakk,
						TMAX, KK, FFC);

				ArrayList<Double> degs = new ArrayList<Double>();// for
																	// calculate
																	// the
																	// degre
																	// of
																	// each
																	// mass
																	// released

				double Mran, MCran, TRran, Xran, openran, outran;

				Mran = mon.simule(TrianM.getA(), TrianM.getC(), monteCarlo);

				if (TrianM.getA() != TrianM.getC())
					degs.add(TrianM.getMembershipDegree(Mran));

				MCran = mon.simule(TrianMC.getA(), TrianMC.getC(), monteCarlo);

				if (TrianMC.getA() != TrianMC.getC())
					degs.add(TrianMC.getMembershipDegree(MCran));

				TRran = mon.simule(TrianTR.getA(), TrianTR.getC(), monteCarlo);

				if (TrianTR.getA() != TrianTR.getC())
					degs.add(TrianTR.getMembershipDegree(TRran));

				Xran = mon.simule(TrianX.getA(), TrianX.getC(), monteCarlo);
				if (TrianX.getA() != TrianX.getC())
					degs.add(TrianX.getMembershipDegree(Xran));

				openran = mon.simule(TrianOpen.getA(), TrianOpen.getC(),
						monteCarlo);
				if (TrianOpen.getA() != TrianOpen.getC())
					degs.add(TrianOpen.getMembershipDegree(openran));

				outran = mon.simule(TrianOut.getA(), TrianOut.getC(),
						monteCarlo);
				if (TrianOut.getA() != TrianOut.getC())
					degs.add(TrianOut.getMembershipDegree(outran));// the random
																	// variables
																	// of
																	// each
																	// tuple

				// sort the array list
				Collections.sort(degs);

				double[] y = new double[] { Mran, MCran, TRran, Xran, openran,
						outran }; // initial state
				double[] ydot = new double[6];

				StepHandler stepHandler = new StepHandler() {

					public void init(double t0, double[] y0, double t) {

					}

					public void handleStep(StepInterpolator interpolator,
							boolean isLast) {
						double t = interpolator.getCurrentTime();
						double[] y = interpolator.getInterpolatedState();

						somme += y[5];

					}
				};

				dp853.addStepHandler(stepHandler);
				dp853.integrate(ode, 0.0, y, 800, y);
				System.out.println(counter + 1);
				dp853.clearEventHandlers();

				MassRelacheeList.add(somme);
				DegreeOfMass.add(degs.get(0));
				MassDegre.put(somme, degs.get(0));
				DegreeMass.put(degs.get(0), somme);
				lblMassReleased1.setText("    " + lblMassReleased1.getText()
						+ "\n" + somme + "  " + degs.get(0));
			}
			Collections.sort(MassRelacheeList);// for calculate the max of
												// mass relased
			MAXmassReleased = MassRelacheeList.get(MassRelacheeList.size() - 1);
			MINmassReleased = MassRelacheeList.get(0);
			MAXdeg = DegreeOfMass.get(DegreeOfMass.size() - 1);
			MOYmassReleased = DegreeMass.get(MAXdeg);
			lblMaxMass.setText("" + MAXmassReleased);
		}

	}

	/**
	 * Create the frame.
	 */
	public FuzzyParamaterWithFixedTimeFailureMain() {

		groupemon = new ButtonGroup();
		Normal = new JCheckBox("Uniform");
		Normal.setBounds(86, 208, 75, 23);
		groupemon.add(Normal);
		Normal.setSelected(true);

		Gaussian = new JCheckBox("Gaussian");
		Gaussian.setBounds(163, 208, 86, 23);
		groupemon.add(Gaussian);
		MC.setColumns(10);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 868, 500);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);

		lamda_iLabel = new JLabel();
		String[] columnName = { "" };
		Object[][] data = { { "masse relache" }, { "" }, { "rendement" },
				{ "" }, { "Lamda_i" }, { "" }, { "delta_H" } };
		Ren = new ButtonGroup();
		mass = new ButtonGroup();
		JPanel ecran = new JPanel();
		ecran.setBorder(BorderFactory
				.createTitledBorder("simulation of the reactor"));
		JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();
		// ecran.setBackground(Color.cyan);

		JPanel panel = new JPanel();
		panel.setBorder(BorderFactory
				.createTitledBorder("Result of Mass Released"));

		panel_2 = new JPanel();
		panel_2.setBorder(BorderFactory
				.createTitledBorder("Calcul of Distance: "));
		GroupLayout gl_contentPane = new GroupLayout(contentPane);
		gl_contentPane
				.setHorizontalGroup(gl_contentPane
						.createParallelGroup(Alignment.LEADING)
						.addGroup(
								gl_contentPane
										.createSequentialGroup()
										.addContainerGap()
										.addGroup(
												gl_contentPane
														.createParallelGroup(
																Alignment.TRAILING,
																false)
														.addComponent(
																panel_2,
																Alignment.LEADING,
																0, 0,
																Short.MAX_VALUE)
														.addGroup(
																Alignment.LEADING,
																gl_contentPane
																		.createSequentialGroup()
																		.addComponent(
																				ecran,
																				GroupLayout.PREFERRED_SIZE,
																				519,
																				GroupLayout.PREFERRED_SIZE)
																		.addPreferredGap(
																				ComponentPlacement.UNRELATED)
																		.addComponent(
																				panel,
																				GroupLayout.PREFERRED_SIZE,
																				294,
																				GroupLayout.PREFERRED_SIZE)))
										.addContainerGap(205, Short.MAX_VALUE)));
		gl_contentPane
				.setVerticalGroup(gl_contentPane
						.createParallelGroup(Alignment.LEADING)
						.addGroup(
								gl_contentPane
										.createSequentialGroup()
										.addComponent(panel_2,
												GroupLayout.PREFERRED_SIZE,
												126, GroupLayout.PREFERRED_SIZE)
										.addPreferredGap(
												ComponentPlacement.UNRELATED)
										.addGroup(
												gl_contentPane
														.createParallelGroup(
																Alignment.TRAILING)
														.addComponent(
																ecran,
																Alignment.LEADING,
																GroupLayout.PREFERRED_SIZE,
																302,
																GroupLayout.PREFERRED_SIZE)
														.addComponent(
																panel,
																GroupLayout.DEFAULT_SIZE,
																303,
																Short.MAX_VALUE))
										.addContainerGap()));

		JLabel lblNewLabel_1 = new JLabel("Rendement");
		lblNewLabel_1.setHorizontalAlignment(SwingConstants.RIGHT);

		renA = new JTextField();
		renA.setText("5");
		renA.setColumns(10);

		renB = new JTextField();
		renB.setText("7.5");
		renB.setColumns(10);

		renC = new JTextField();
		renC.setText("10");
		renC.setColumns(10);

		FuzzyRen = new JRadioButton("Fuzzy_Number");
		Ren.add(FuzzyRen);

		NormalRen = new JRadioButton("Normal_Number");
		Ren.add(NormalRen);
		FuzzyRen.setSelected(true);
		JLabel lblNewLabel_2 = new JLabel("Lamda_i");
		lblNewLabel_2.setHorizontalAlignment(SwingConstants.RIGHT);

		lamda = new JTextField();
		lamda.setText("22");
		lamda.setColumns(10);

		JLabel lblNewLabel_4 = new JLabel(
				"Distance of dangers using TNTequivalent");
		lblNewLabel_4.setHorizontalAlignment(SwingConstants.LEFT);

		JLabel lblNewLabel_3 = new JLabel("delta_H");
		lblNewLabel_3.setHorizontalAlignment(SwingConstants.RIGHT);

		deltah = new JTextField();
		deltah.setText("48220");
		deltah.setColumns(10);

		JButton calcul = new JButton("click here");
		calcul.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				// TODO Auto-generated method stub
				ResultDis.setText("");
				new DistanceThread().start();
			}
		});
		JLabel lblResult = new JLabel("Result =");
		lblResult.setForeground(Color.RED);
		lblResult.setFont(new Font("Tahoma", Font.BOLD, 11));

		ResultDis = new JLabel("");
		ResultDis.setBackground(Color.RED);
		GroupLayout gl_panel_2 = new GroupLayout(panel_2);
		gl_panel_2
				.setHorizontalGroup(gl_panel_2
						.createParallelGroup(Alignment.LEADING)
						.addGroup(
								gl_panel_2
										.createSequentialGroup()
										.addGap(14)
										.addGroup(
												gl_panel_2
														.createParallelGroup(
																Alignment.LEADING)
														.addGroup(
																gl_panel_2
																		.createSequentialGroup()
																		.addComponent(
																				lblNewLabel_1)
																		.addPreferredGap(
																				ComponentPlacement.RELATED)
																		.addComponent(
																				renA,
																				GroupLayout.PREFERRED_SIZE,
																				GroupLayout.DEFAULT_SIZE,
																				GroupLayout.PREFERRED_SIZE)
																		.addGap(6)
																		.addComponent(
																				renB,
																				GroupLayout.PREFERRED_SIZE,
																				GroupLayout.DEFAULT_SIZE,
																				GroupLayout.PREFERRED_SIZE)
																		.addGap(6)
																		.addComponent(
																				renC,
																				GroupLayout.PREFERRED_SIZE,
																				GroupLayout.DEFAULT_SIZE,
																				GroupLayout.PREFERRED_SIZE)
																		.addGap(6)
																		.addComponent(
																				FuzzyRen,
																				GroupLayout.PREFERRED_SIZE,
																				113,
																				GroupLayout.PREFERRED_SIZE)
																		.addGap(2)
																		.addComponent(
																				NormalRen,
																				GroupLayout.PREFERRED_SIZE,
																				142,
																				GroupLayout.PREFERRED_SIZE))
														.addGroup(
																gl_panel_2
																		.createSequentialGroup()
																		.addGap(16)
																		.addGroup(
																				gl_panel_2
																						.createParallelGroup(
																								Alignment.TRAILING)
																						.addComponent(
																								lblNewLabel_2)
																						.addComponent(
																								lblNewLabel_3))
																		.addPreferredGap(
																				ComponentPlacement.RELATED)
																		.addGroup(
																				gl_panel_2
																						.createParallelGroup(
																								Alignment.LEADING)
																						.addGroup(
																								gl_panel_2
																										.createSequentialGroup()
																										.addComponent(
																												lamda,
																												GroupLayout.PREFERRED_SIZE,
																												GroupLayout.DEFAULT_SIZE,
																												GroupLayout.PREFERRED_SIZE)
																										.addPreferredGap(
																												ComponentPlacement.UNRELATED)
																										.addComponent(
																												lblNewLabel_4,
																												GroupLayout.PREFERRED_SIZE,
																												254,
																												GroupLayout.PREFERRED_SIZE))
																						.addGroup(
																								gl_panel_2
																										.createSequentialGroup()
																										.addComponent(
																												deltah,
																												GroupLayout.PREFERRED_SIZE,
																												GroupLayout.DEFAULT_SIZE,
																												GroupLayout.PREFERRED_SIZE)
																										.addGap(8)
																										.addComponent(
																												calcul,
																												GroupLayout.PREFERRED_SIZE,
																												97,
																												GroupLayout.PREFERRED_SIZE)
																										.addPreferredGap(
																												ComponentPlacement.RELATED)
																										.addComponent(
																												lblResult)
																										.addPreferredGap(
																												ComponentPlacement.RELATED)
																										.addComponent(
																												ResultDis,
																												GroupLayout.PREFERRED_SIZE,
																												430,
																												GroupLayout.PREFERRED_SIZE)))))
										.addContainerGap(117, Short.MAX_VALUE)));
		gl_panel_2
				.setVerticalGroup(gl_panel_2
						.createParallelGroup(Alignment.LEADING)
						.addGroup(
								gl_panel_2
										.createSequentialGroup()
										.addContainerGap()
										.addGroup(
												gl_panel_2
														.createParallelGroup(
																Alignment.LEADING)
														.addGroup(
																gl_panel_2
																		.createSequentialGroup()
																		.addGap(12)
																		.addGroup(
																				gl_panel_2
																						.createParallelGroup(
																								Alignment.BASELINE)
																						.addComponent(
																								renA,
																								GroupLayout.PREFERRED_SIZE,
																								GroupLayout.DEFAULT_SIZE,
																								GroupLayout.PREFERRED_SIZE)
																						.addComponent(
																								lblNewLabel_1)))
														.addGroup(
																gl_panel_2
																		.createSequentialGroup()
																		.addGap(12)
																		.addComponent(
																				renB,
																				GroupLayout.PREFERRED_SIZE,
																				GroupLayout.DEFAULT_SIZE,
																				GroupLayout.PREFERRED_SIZE))
														.addGroup(
																gl_panel_2
																		.createSequentialGroup()
																		.addGap(12)
																		.addComponent(
																				renC,
																				GroupLayout.PREFERRED_SIZE,
																				GroupLayout.DEFAULT_SIZE,
																				GroupLayout.PREFERRED_SIZE))
														.addGroup(
																gl_panel_2
																		.createSequentialGroup()
																		.addGap(11)
																		.addComponent(
																				FuzzyRen))
														.addGroup(
																gl_panel_2
																		.createSequentialGroup()
																		.addGap(11)
																		.addComponent(
																				NormalRen)))
										.addPreferredGap(
												ComponentPlacement.UNRELATED)
										.addGroup(
												gl_panel_2
														.createParallelGroup(
																Alignment.LEADING)
														.addGroup(
																gl_panel_2
																		.createParallelGroup(
																				Alignment.BASELINE)
																		.addComponent(
																				lamda,
																				GroupLayout.PREFERRED_SIZE,
																				GroupLayout.DEFAULT_SIZE,
																				GroupLayout.PREFERRED_SIZE)
																		.addComponent(
																				lblNewLabel_2))
														.addComponent(
																lblNewLabel_4,
																GroupLayout.PREFERRED_SIZE,
																19,
																GroupLayout.PREFERRED_SIZE))
										.addPreferredGap(
												ComponentPlacement.RELATED)
										.addGroup(
												gl_panel_2
														.createParallelGroup(
																Alignment.BASELINE)
														.addComponent(calcul)
														.addComponent(lblResult)
														.addComponent(
																deltah,
																GroupLayout.PREFERRED_SIZE,
																GroupLayout.DEFAULT_SIZE,
																GroupLayout.PREFERRED_SIZE)
														.addComponent(
																lblNewLabel_3)
														.addComponent(ResultDis))
										.addContainerGap(34, Short.MAX_VALUE)));
		panel_2.setLayout(gl_panel_2);

		lblCounter = new JLabel("");
		lblCounter_1 = new JLabel("Counter");
		lblMassReleased1 = new JTextArea();
		scr = new JScrollPane(lblMassReleased1);
		panel.add(scr);
		scr.setBounds(15, 55, 270, 200);

		JLabel lblNewLabel_6 = new JLabel("MAXmass_released =");

		lblMaxMass = new JLabel("");
		GroupLayout gl_panel = new GroupLayout(panel);
		gl_panel.setHorizontalGroup(gl_panel
				.createParallelGroup(Alignment.TRAILING)
				.addGroup(
						gl_panel.createSequentialGroup()
								.addContainerGap()
								.addGroup(
										gl_panel.createParallelGroup(
												Alignment.LEADING)
												.addGroup(
														gl_panel.createSequentialGroup()
																.addComponent(
																		lblCounter_1)
																.addPreferredGap(
																		ComponentPlacement.UNRELATED)
																.addComponent(
																		lblCounter))
												.addGroup(
														gl_panel.createSequentialGroup()
																.addComponent(
																		lblNewLabel_6)
																.addPreferredGap(
																		ComponentPlacement.RELATED)
																.addComponent(
																		lblMaxMass,
																		GroupLayout.DEFAULT_SIZE,
																		144,
																		Short.MAX_VALUE)))
								.addContainerGap()));
		gl_panel.setVerticalGroup(gl_panel.createParallelGroup(
				Alignment.LEADING)
				.addGroup(
						gl_panel.createSequentialGroup()
								.addContainerGap()
								.addGroup(
										gl_panel.createParallelGroup(
												Alignment.BASELINE)
												.addComponent(lblCounter_1)
												.addComponent(lblCounter))
								.addPreferredGap(ComponentPlacement.RELATED,
										230, Short.MAX_VALUE)
								.addGroup(
										gl_panel.createParallelGroup(
												Alignment.BASELINE)
												.addComponent(lblNewLabel_6)
												.addComponent(lblMaxMass))
								.addContainerGap()));
		panel.setLayout(gl_panel);
		ButtonGroup groupemon = new ButtonGroup();

		panel_3 = new JPanel();
		panel_3.setBorder(BorderFactory.createTitledBorder("Initial State"));
		GroupLayout gl_ecran = new GroupLayout(ecran);
		gl_ecran.setHorizontalGroup(gl_ecran.createParallelGroup(
				Alignment.LEADING).addComponent(panel_3,
				GroupLayout.DEFAULT_SIZE, 507, Short.MAX_VALUE));
		gl_ecran.setVerticalGroup(gl_ecran.createParallelGroup(
				Alignment.LEADING).addGroup(
				gl_ecran.createSequentialGroup()
						.addContainerGap()
						.addComponent(panel_3, GroupLayout.PREFERRED_SIZE, 252,
								GroupLayout.PREFERRED_SIZE)
						.addContainerGap(16, Short.MAX_VALUE)));

		JLabel lblTr = new JLabel("TR");
		lblTr.setBounds(31, 92, 45, 14);
		lblTr.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblX = new JLabel("X");
		lblX.setBounds(31, 118, 45, 14);
		lblX.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblMc = new JLabel("MC");
		lblMc.setBounds(31, 66, 45, 14);
		lblMc.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblM = new JLabel("M");
		lblM.setBounds(31, 37, 45, 14);
		lblM.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblOpen = new JLabel("Open");
		lblOpen.setBounds(31, 144, 45, 14);
		lblOpen.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblOut = new JLabel("out");
		lblOut.setBounds(31, 167, 45, 14);
		lblOut.setHorizontalAlignment(SwingConstants.RIGHT);

		out = new JTextField();
		out.setBounds(86, 164, 45, 20);
		out.setText("0");
		out.setColumns(10);

		open = new JTextField();
		open.setText("0");
		open.setBounds(86, 141, 45, 20);
		open.setColumns(10);

		X = new JTextField();
		X.setText("0");
		X.setBounds(86, 115, 45, 20);
		X.setColumns(10);

		TR = new JTextField();
		TR.setText("80");
		TR.setBounds(86, 89, 45, 20);
		TR.setColumns(10);

		M = new JTextField();
		M.setText("4400");
		M.setBounds(86, 37, 45, 20);
		M.setColumns(10);

		Mc = new JTextField();
		Mc.setText("0");
		Mc.setBounds(86, 63, 45, 20);
		Mc.setColumns(10);

		M1 = new JTextField();
		M1.setText("5000");
		M1.setBounds(137, 37, 45, 20);
		M1.setColumns(10);

		out1 = new JTextField();
		out1.setText("0");
		out1.setBounds(137, 164, 45, 20);
		out1.setColumns(10);

		open1 = new JTextField();
		open1.setText("0");
		open1.setBounds(137, 141, 45, 20);
		open1.setColumns(10);

		X1 = new JTextField();
		X1.setText("0");
		X1.setBounds(137, 115, 45, 20);
		X1.setColumns(10);

		TR1 = new JTextField();
		TR1.setText("90");
		TR1.setBounds(137, 89, 45, 20);
		TR1.setColumns(10);

		Mc1 = new JTextField();
		Mc1.setText("0");
		Mc1.setBounds(137, 63, 45, 20);
		Mc1.setColumns(10);

		M2 = new JTextField();
		M2.setText("6000");
		M2.setBounds(188, 37, 45, 20);
		M2.setColumns(10);

		Mc2 = new JTextField();
		Mc2.setText("0");
		Mc2.setBounds(188, 63, 45, 20);
		Mc2.setColumns(10);

		TR2 = new JTextField();
		TR2.setText("100");
		TR2.setBounds(188, 89, 45, 20);
		TR2.setColumns(10);

		X2 = new JTextField();
		X2.setText("0");
		X2.setBounds(188, 115, 45, 20);
		X2.setColumns(10);

		open2 = new JTextField();
		open2.setText("0");
		open2.setBounds(188, 141, 45, 20);
		open2.setColumns(10);

		out2 = new JTextField();
		out2.setText("0");
		out2.setBounds(188, 164, 45, 20);
		out2.setColumns(10);

		JButton btnNewButton = new JButton("Simulation");
		btnNewButton.setBounds(257, 208, 102, 23);
		btnNewButton.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				// TODO Auto-generated method stub
				new MyThread().start();
			}
		});

		JLabel lblTbreakdown = new JLabel("Tbreakdown");
		lblTbreakdown.setBounds(293, 40, 86, 14);
		lblTbreakdown.setHorizontalAlignment(SwingConstants.RIGHT);

		tbreak = new JTextField();
		tbreak.setText("700");
		tbreak.setBounds(383, 37, 86, 20);
		tbreak.setColumns(10);

		JLabel lblNewLabel_5 = new JLabel("Tmax");
		lblNewLabel_5.setBounds(293, 66, 86, 14);
		lblNewLabel_5.setHorizontalAlignment(SwingConstants.RIGHT);

		tmax = new JTextField();
		tmax.setBounds(383, 63, 86, 20);
		tmax.setText("712");
		tmax.setColumns(10);

		JLabel lblK = new JLabel("variation of k");
		lblK.setBounds(293, 92, 86, 14);
		lblK.setHorizontalAlignment(SwingConstants.RIGHT);

		k = new JTextField();
		k.setBounds(383, 89, 86, 20);
		k.setText("1");
		k.setColumns(10);

		JLabel lblFc = new JLabel("Fc");
		lblFc.setBounds(293, 118, 86, 14);
		lblFc.setHorizontalAlignment(SwingConstants.RIGHT);

		FcValue = new JTextField();
		FcValue.setText("3300");
		FcValue.setBounds(383, 115, 86, 20);
		FcValue.setColumns(10);

		lblLoi = new JLabel("Loi");
		lblLoi.setBounds(31, 212, 45, 14);
		lblLoi.setHorizontalAlignment(SwingConstants.RIGHT);
		panel_3.setLayout(null);
		panel_3.add(lblM);
		panel_3.add(M);
		panel_3.add(M1);
		panel_3.add(M2);
		panel_3.add(lblTbreakdown);
		panel_3.add(tbreak);
		panel_3.add(lblMc);
		panel_3.add(Mc);
		panel_3.add(Mc1);
		panel_3.add(Mc2);
		panel_3.add(lblNewLabel_5);
		panel_3.add(tmax);
		panel_3.add(lblTr);
		panel_3.add(TR);
		panel_3.add(TR1);
		panel_3.add(TR2);
		panel_3.add(lblK);
		panel_3.add(k);
		panel_3.add(lblX);
		panel_3.add(X);
		panel_3.add(X1);
		panel_3.add(X2);
		panel_3.add(lblFc);
		panel_3.add(FcValue);
		panel_3.add(lblOpen);
		panel_3.add(open);
		panel_3.add(open1);
		panel_3.add(open2);
		panel_3.add(btnNewButton);
		panel_3.add(lblOut);
		panel_3.add(out);
		panel_3.add(out1);
		panel_3.add(out2);
		panel_3.add(lblLoi);
		panel_3.add(Normal);
		panel_3.add(Gaussian);

		JLabel lblNewLabel = new JLabel("Sample Nb.");
		lblNewLabel.setHorizontalAlignment(SwingConstants.RIGHT);
		lblNewLabel.setBounds(293, 167, 86, 14);
		panel_3.add(lblNewLabel);

		txtlSampleNb = new JTextField();
		txtlSampleNb.setText("1000");
		txtlSampleNb.setBounds(383, 164, 86, 20);
		panel_3.add(txtlSampleNb);
		txtlSampleNb.setColumns(10);

		textCut = new JTextField();
		textCut.setText("1");
		textCut.setBounds(383, 141, 86, 20);
		panel_3.add(textCut);
		textCut.setColumns(10);

		Cut = new JLabel("Cut");
		Cut.setHorizontalAlignment(SwingConstants.RIGHT);
		Cut.setBounds(332, 143, 46, 14);
		panel_3.add(Cut);
		ecran.setLayout(gl_ecran);
		gl_contentPane.setHonorsVisibility(false);
		contentPane.setLayout(gl_contentPane);
	}
}
