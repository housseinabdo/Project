import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Image;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.plaf.basic.BasicBorders.RadioButtonBorder;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.GroupLayout;
import javax.swing.GroupLayout.Alignment;
import javax.swing.JTextField;

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

import classes.MonteCarlo;
import classes.NormalInputODE;

import com.panayotis.gnuplot.JavaPlot;
import com.panayotis.gnuplot.plot.AbstractPlot;
import com.panayotis.gnuplot.plot.DataSetPlot;
import com.panayotis.gnuplot.style.PlotStyle;
import com.panayotis.gnuplot.style.Style;

import java.awt.Font;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import classes.*;
public class NormalinputMain extends JFrame {

	private JPanel contentPane;
	public JTextField mA;
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
	private JLabel massReleased;
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

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {

		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					NormalinputMain frame = new NormalinputMain();
					frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	public NormalinputMain() {
		MC.setColumns(10);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(100, 100, 800, 500);
		contentPane = new JPanel();
		contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
		setContentPane(contentPane);

		mA = new JTextField();
		mA.setColumns(10);

		lamda_iLabel = new JLabel();

		renA = new JTextField();
		renA.setColumns(10);

		lamda = new JTextField();
		lamda.setColumns(10);

		renB = new JTextField();
		renB.setColumns(10);

		deltah = new JTextField();
		deltah.setColumns(10);

		JButton calcul = new JButton("click here");
		calcul.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent arg0) {
				// TODO Auto-generated method stub
				double lamdai, m, deltaH;
				lamdai = Double.parseDouble(lamda.getText());
				m = Double.parseDouble(mA.getText());
				deltaH = Double.parseDouble(deltah.getText());
				if (FuzzyRen.isSelected()) {
					try {
						ResultDis.setText("");
						TriangularFuzzyNumber rendement = new TriangularFuzzyNumber(
								Double.parseDouble(renA.getText()), Double
										.parseDouble(renB.getText()), Double
										.parseDouble(renC.getText()));

						double a = lamdai
								* (rendement.getA() * ((m * deltaH) / 4690));
						double b = lamdai
								* (rendement.getB() * ((m * deltaH) / 4690));
						double c = lamdai
								* (rendement.getC() * ((m * deltaH) / 4690));
						TriangularFuzzyNumber distance = new TriangularFuzzyNumber(
								a, b, c);
						// System.out.println(distance);

						ResultDis.setText("(" + distance.getA() + " , "
								+ distance.getB() + " , " + distance.getC()
								+ ")");

					} catch (FuzzyOntologyException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				if (NormalRen.isSelected()) {
					double distan = lamdai
							* (Double.parseDouble(renA.getText()) * ((m * deltaH) / 4690));
					// System.out.println(distan);
					ResultDis.setText("distance = " + distan);
				}
				mA.setText("");
				renA.setText("");
				renB.setText("");
				renC.setText("");
				lamda.setText("");
				deltah.setText("");
			}
		});
		String[] columnName = { "" };
		Object[][] data = { { "masse relache" }, { "" }, { "rendement" },
				{ "" }, { "Lamda_i" }, { "" }, { "delta_H" } };

		renC = new JTextField();
		renC.setColumns(10);

		JLabel lblNewLabel = new JLabel("mass releases");
		lblNewLabel.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblNewLabel_1 = new JLabel("Rendement");
		lblNewLabel_1.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblNewLabel_2 = new JLabel("Lamda_i");
		lblNewLabel_2.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblNewLabel_3 = new JLabel("delta_H");
		lblNewLabel_3.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblNewLabel_4 = new JLabel(
				"Distance of dangers using TNTequivalent");
		lblNewLabel_4.setHorizontalAlignment(SwingConstants.LEFT);

		FuzzyRen = new JRadioButton("Fuzzy_Number");

		NormalRen = new JRadioButton("Normal_Number");
		Ren = new ButtonGroup();
		mass = new ButtonGroup();
		Ren.add(FuzzyRen);
		Ren.add(NormalRen);
		JPanel ecran = new JPanel();

		ecran.setBackground(Color.cyan);
		JLabel lblResult = new JLabel("Result =");
		lblResult.setForeground(Color.RED);
		lblResult.setFont(new Font("Tahoma", Font.BOLD, 11));

		ResultDis = new JLabel("");

		JPanel panel = new JPanel();
		panel.setBackground(Color.black);
		GroupLayout gl_contentPane = new GroupLayout(contentPane);
		gl_contentPane
				.setHorizontalGroup(gl_contentPane
						.createParallelGroup(Alignment.LEADING)
						.addGroup(
								gl_contentPane
										.createSequentialGroup()
										.addGroup(
												gl_contentPane
														.createParallelGroup(
																Alignment.LEADING)
														.addGroup(
																gl_contentPane
																		.createSequentialGroup()
																		.addGap(20)
																		.addGroup(
																				gl_contentPane
																						.createParallelGroup(
																								Alignment.LEADING,
																								false)
																						.addComponent(
																								lblNewLabel,
																								GroupLayout.DEFAULT_SIZE,
																								87,
																								Short.MAX_VALUE)
																						.addComponent(
																								lblNewLabel_1,
																								GroupLayout.DEFAULT_SIZE,
																								GroupLayout.DEFAULT_SIZE,
																								Short.MAX_VALUE)
																						.addComponent(
																								lblNewLabel_2,
																								GroupLayout.DEFAULT_SIZE,
																								GroupLayout.DEFAULT_SIZE,
																								Short.MAX_VALUE)
																						.addComponent(
																								lblNewLabel_3,
																								GroupLayout.DEFAULT_SIZE,
																								GroupLayout.DEFAULT_SIZE,
																								Short.MAX_VALUE))
																		.addGap(18)
																		.addGroup(
																				gl_contentPane
																						.createParallelGroup(
																								Alignment.LEADING)
																						.addComponent(
																								mA,
																								GroupLayout.PREFERRED_SIZE,
																								GroupLayout.DEFAULT_SIZE,
																								GroupLayout.PREFERRED_SIZE)
																						.addComponent(
																								renA,
																								GroupLayout.PREFERRED_SIZE,
																								GroupLayout.DEFAULT_SIZE,
																								GroupLayout.PREFERRED_SIZE)
																						.addComponent(
																								lamda,
																								GroupLayout.PREFERRED_SIZE,
																								GroupLayout.DEFAULT_SIZE,
																								GroupLayout.PREFERRED_SIZE)
																						.addComponent(
																								deltah,
																								GroupLayout.PREFERRED_SIZE,
																								GroupLayout.DEFAULT_SIZE,
																								GroupLayout.PREFERRED_SIZE))
																		.addGap(6)
																		.addGroup(
																				gl_contentPane
																						.createParallelGroup(
																								Alignment.LEADING)
																						.addComponent(
																								lblNewLabel_4,
																								GroupLayout.PREFERRED_SIZE,
																								260,
																								GroupLayout.PREFERRED_SIZE)
																						.addGroup(
																								gl_contentPane
																										.createSequentialGroup()
																										.addComponent(
																												renB,
																												GroupLayout.PREFERRED_SIZE,
																												GroupLayout.DEFAULT_SIZE,
																												GroupLayout.PREFERRED_SIZE)
																										.addPreferredGap(
																												ComponentPlacement.RELATED)
																										.addComponent(
																												renC,
																												GroupLayout.PREFERRED_SIZE,
																												GroupLayout.DEFAULT_SIZE,
																												GroupLayout.PREFERRED_SIZE)
																										.addPreferredGap(
																												ComponentPlacement.UNRELATED)
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
																								gl_contentPane
																										.createSequentialGroup()
																										.addComponent(
																												calcul,
																												GroupLayout.PREFERRED_SIZE,
																												97,
																												GroupLayout.PREFERRED_SIZE)
																										.addPreferredGap(
																												ComponentPlacement.UNRELATED)
																										.addComponent(
																												lblResult)
																										.addPreferredGap(
																												ComponentPlacement.RELATED)
																										.addComponent(
																												ResultDis))))
														.addGroup(
																gl_contentPane
																		.createSequentialGroup()
																		.addGap(35)
																		.addComponent(
																				ecran,
																				GroupLayout.PREFERRED_SIZE,
																				GroupLayout.DEFAULT_SIZE,
																				GroupLayout.PREFERRED_SIZE)
																		.addGap(18)
																		.addComponent(
																				panel,
																				GroupLayout.DEFAULT_SIZE,
																				291,
																				Short.MAX_VALUE)))
										.addContainerGap()));
		gl_contentPane
				.setVerticalGroup(gl_contentPane
						.createParallelGroup(Alignment.LEADING)
						.addGroup(
								gl_contentPane
										.createSequentialGroup()
										.addGroup(
												gl_contentPane
														.createParallelGroup(
																Alignment.BASELINE)
														.addComponent(
																mA,
																GroupLayout.PREFERRED_SIZE,
																GroupLayout.DEFAULT_SIZE,
																GroupLayout.PREFERRED_SIZE)
														.addComponent(
																lblNewLabel))
										.addGap(11)
										.addGroup(
												gl_contentPane
														.createParallelGroup(
																Alignment.BASELINE)
														.addComponent(
																renA,
																GroupLayout.PREFERRED_SIZE,
																GroupLayout.DEFAULT_SIZE,
																GroupLayout.PREFERRED_SIZE)
														.addComponent(
																lblNewLabel_1)
														.addComponent(
																renB,
																GroupLayout.PREFERRED_SIZE,
																GroupLayout.DEFAULT_SIZE,
																GroupLayout.PREFERRED_SIZE)
														.addComponent(
																renC,
																GroupLayout.PREFERRED_SIZE,
																GroupLayout.DEFAULT_SIZE,
																GroupLayout.PREFERRED_SIZE)
														.addComponent(FuzzyRen)
														.addComponent(NormalRen))
										.addGap(6)
										.addGroup(
												gl_contentPane
														.createParallelGroup(
																Alignment.LEADING)
														.addGroup(
																gl_contentPane
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
										.addGap(11)
										.addGroup(
												gl_contentPane
														.createParallelGroup(
																Alignment.BASELINE)
														.addComponent(
																deltah,
																GroupLayout.PREFERRED_SIZE,
																GroupLayout.DEFAULT_SIZE,
																GroupLayout.PREFERRED_SIZE)
														.addComponent(
																lblNewLabel_3)
														.addComponent(calcul)
														.addComponent(lblResult)
														.addComponent(ResultDis))
										.addPreferredGap(
												ComponentPlacement.RELATED, 18,
												Short.MAX_VALUE)
										.addGroup(
												gl_contentPane
														.createParallelGroup(
																Alignment.LEADING,
																false)
														.addComponent(
																panel,
																GroupLayout.DEFAULT_SIZE,
																GroupLayout.DEFAULT_SIZE,
																Short.MAX_VALUE)
														.addComponent(
																ecran,
																GroupLayout.DEFAULT_SIZE,
																300,
																Short.MAX_VALUE))
										.addGap(19)));

		massReleased = new JLabel("");
		massReleased.setForeground(Color.RED);
		massReleased.setFont(new Font("Tahoma", Font.PLAIN, 30));
		GroupLayout gl_panel = new GroupLayout(panel);
		gl_panel.setHorizontalGroup(gl_panel.createParallelGroup(
				Alignment.LEADING).addGroup(
				gl_panel.createSequentialGroup()
						.addContainerGap()
						.addComponent(massReleased, GroupLayout.DEFAULT_SIZE,
								271, Short.MAX_VALUE).addContainerGap()));
		gl_panel.setVerticalGroup(gl_panel.createParallelGroup(
				Alignment.LEADING).addGroup(
				Alignment.TRAILING,
				gl_panel.createSequentialGroup()
						.addContainerGap()
						.addComponent(massReleased, GroupLayout.DEFAULT_SIZE,
								278, Short.MAX_VALUE).addContainerGap()));
		panel.setLayout(gl_panel);

		JButton btnNewButton = new JButton("Simulation");

		JLabel lblSimulationOfThe = new JLabel("simulation of the reactor");
		lblSimulationOfThe.setFont(new Font("Tahoma", Font.PLAIN, 15));
		lblSimulationOfThe.setForeground(Color.BLUE);

		JLabel lblInitialStates = new JLabel("Initial states");
		lblInitialStates.setHorizontalAlignment(SwingConstants.CENTER);

		JLabel lblM = new JLabel("M");
		lblM.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblMc = new JLabel("MC");
		lblMc.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblTr = new JLabel("TR");
		lblTr.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblX = new JLabel("X");
		lblX.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblOpen = new JLabel("Open");
		lblOpen.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblOut = new JLabel("out");
		lblOut.setHorizontalAlignment(SwingConstants.RIGHT);

		M = new JTextField();
		M.setText("4400");
		M.setColumns(10);

		Mc = new JTextField();
		Mc.setText("0");
		Mc.setColumns(10);

		TR = new JTextField();
		TR.setText("80");
		TR.setColumns(10);

		X = new JTextField();
		X.setText("0");
		X.setColumns(10);

		open = new JTextField();
		open.setText("0");
		open.setColumns(10);

		out = new JTextField();
		out.setText("0");
		out.setColumns(10);

		JLabel lblTbreakdown = new JLabel("Tbreakdown");
		lblTbreakdown.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblNewLabel_5 = new JLabel("Tmax");
		lblNewLabel_5.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblK = new JLabel("variation of k");
		lblK.setHorizontalAlignment(SwingConstants.RIGHT);

		JLabel lblFc = new JLabel("Fc");
		lblFc.setHorizontalAlignment(SwingConstants.RIGHT);

		tbreak = new JTextField();
		tbreak.setText("700");
		tbreak.setColumns(10);

		tmax = new JTextField();
		tmax.setText("");
		tmax.setColumns(10);

		k = new JTextField();
		k.setText("1");
		k.setColumns(10);

		FcValue = new JTextField();
		FcValue.setText("3300");
		FcValue.setColumns(10);
		GroupLayout gl_ecran = new GroupLayout(ecran);
		gl_ecran.setHorizontalGroup(gl_ecran
				.createParallelGroup(Alignment.TRAILING)
				.addGroup(
						gl_ecran.createSequentialGroup().addGap(139)
								.addComponent(lblSimulationOfThe)
								.addContainerGap(123, Short.MAX_VALUE))
				.addGroup(
						gl_ecran.createSequentialGroup()
								.addGroup(
										gl_ecran.createParallelGroup(
												Alignment.LEADING)
												.addGroup(
														gl_ecran.createSequentialGroup()
																.addContainerGap()
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.TRAILING)
																				.addComponent(
																						lblOut,
																						Alignment.LEADING,
																						GroupLayout.DEFAULT_SIZE,
																						26,
																						Short.MAX_VALUE)
																				.addComponent(
																						lblOpen,
																						Alignment.LEADING,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						Short.MAX_VALUE)
																				.addComponent(
																						lblX,
																						Alignment.LEADING,
																						GroupLayout.DEFAULT_SIZE,
																						26,
																						Short.MAX_VALUE)
																				.addComponent(
																						lblTr)
																				.addComponent(
																						lblMc,
																						GroupLayout.DEFAULT_SIZE,
																						26,
																						Short.MAX_VALUE)))
												.addGroup(
														gl_ecran.createSequentialGroup()
																.addGap(18)
																.addComponent(
																		lblM,
																		GroupLayout.DEFAULT_SIZE,
																		18,
																		Short.MAX_VALUE)))
								.addPreferredGap(ComponentPlacement.UNRELATED)
								.addGroup(
										gl_ecran.createParallelGroup(
												Alignment.LEADING)
												.addComponent(
														M,
														GroupLayout.PREFERRED_SIZE,
														GroupLayout.DEFAULT_SIZE,
														GroupLayout.PREFERRED_SIZE)
												.addGroup(
														gl_ecran.createParallelGroup(
																Alignment.TRAILING,
																false)
																.addComponent(
																		btnNewButton,
																		Alignment.LEADING,
																		GroupLayout.DEFAULT_SIZE,
																		GroupLayout.DEFAULT_SIZE,
																		Short.MAX_VALUE)
																.addComponent(
																		out,
																		Alignment.LEADING))
												.addGroup(
														gl_ecran.createSequentialGroup()
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.LEADING)
																				.addGroup(
																						gl_ecran.createParallelGroup(
																								Alignment.LEADING)
																								.addComponent(
																										Mc,
																										GroupLayout.PREFERRED_SIZE,
																										GroupLayout.DEFAULT_SIZE,
																										GroupLayout.PREFERRED_SIZE)
																								.addComponent(
																										TR,
																										GroupLayout.PREFERRED_SIZE,
																										GroupLayout.DEFAULT_SIZE,
																										GroupLayout.PREFERRED_SIZE)
																								.addComponent(
																										X,
																										GroupLayout.PREFERRED_SIZE,
																										GroupLayout.DEFAULT_SIZE,
																										GroupLayout.PREFERRED_SIZE))
																				.addComponent(
																						open,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE))
																.addGap(57)
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.LEADING,
																				false)
																				.addComponent(
																						lblTbreakdown,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						Short.MAX_VALUE)
																				.addComponent(
																						lblNewLabel_5,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						Short.MAX_VALUE)
																				.addComponent(
																						lblK,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						Short.MAX_VALUE)
																				.addComponent(
																						lblFc,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						Short.MAX_VALUE))
																.addPreferredGap(
																		ComponentPlacement.RELATED)
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.LEADING)
																				.addComponent(
																						FcValue,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE)
																				.addComponent(
																						k,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE)
																				.addComponent(
																						tmax,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE)
																				.addComponent(
																						tbreak,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE))))
								.addGap(41))
				.addGroup(
						gl_ecran.createSequentialGroup().addGap(60)
								.addComponent(lblInitialStates)
								.addContainerGap(273, Short.MAX_VALUE)));
		gl_ecran.setVerticalGroup(gl_ecran
				.createParallelGroup(Alignment.LEADING)
				.addGroup(
						gl_ecran.createSequentialGroup()
								.addContainerGap()
								.addGroup(
										gl_ecran.createParallelGroup(
												Alignment.TRAILING)
												.addGroup(
														gl_ecran.createSequentialGroup()
																.addComponent(
																		lblSimulationOfThe)
																.addGap(15)
																.addComponent(
																		lblInitialStates)
																.addPreferredGap(
																		ComponentPlacement.UNRELATED)
																.addComponent(
																		M,
																		GroupLayout.PREFERRED_SIZE,
																		GroupLayout.DEFAULT_SIZE,
																		GroupLayout.PREFERRED_SIZE))
												.addComponent(lblM))
								.addPreferredGap(ComponentPlacement.RELATED)
								.addGroup(
										gl_ecran.createParallelGroup(
												Alignment.BASELINE)
												.addComponent(
														Mc,
														GroupLayout.PREFERRED_SIZE,
														GroupLayout.DEFAULT_SIZE,
														GroupLayout.PREFERRED_SIZE)
												.addComponent(lblMc)
												.addComponent(lblTbreakdown)
												.addComponent(
														tbreak,
														GroupLayout.PREFERRED_SIZE,
														GroupLayout.DEFAULT_SIZE,
														GroupLayout.PREFERRED_SIZE))
								.addPreferredGap(ComponentPlacement.RELATED)
								.addGroup(
										gl_ecran.createParallelGroup(
												Alignment.LEADING)
												.addGroup(
														gl_ecran.createSequentialGroup()
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.BASELINE)
																				.addComponent(
																						TR,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE)
																				.addComponent(
																						lblTr))
																.addPreferredGap(
																		ComponentPlacement.RELATED)
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.BASELINE)
																				.addComponent(
																						X,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE)
																				.addComponent(
																						lblX))
																.addPreferredGap(
																		ComponentPlacement.RELATED)
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.BASELINE)
																				.addComponent(
																						open,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE)
																				.addComponent(
																						lblOpen))
																.addGap(8)
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.BASELINE)
																				.addComponent(
																						out,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE)
																				.addComponent(
																						lblOut))
																.addGap(30)
																.addComponent(
																		btnNewButton))
												.addGroup(
														gl_ecran.createSequentialGroup()
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.BASELINE)
																				.addComponent(
																						lblNewLabel_5)
																				.addComponent(
																						tmax,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE))
																.addPreferredGap(
																		ComponentPlacement.RELATED)
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.BASELINE)
																				.addComponent(
																						lblK)
																				.addComponent(
																						k,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE))
																.addPreferredGap(
																		ComponentPlacement.RELATED)
																.addGroup(
																		gl_ecran.createParallelGroup(
																				Alignment.BASELINE)
																				.addComponent(
																						lblFc)
																				.addComponent(
																						FcValue,
																						GroupLayout.PREFERRED_SIZE,
																						GroupLayout.DEFAULT_SIZE,
																						GroupLayout.PREFERRED_SIZE))))
								.addContainerGap(38, Short.MAX_VALUE)));
		ecran.setLayout(gl_ecran);
		btnNewButton.addActionListener(new ActionListener() {
			double somme = 0;

			@Override
			public void actionPerformed(ActionEvent arg0) {
				// TODO Auto-generated method stub
				somme = 0;
				TriangularFuzzyNumber DEC;
				double t_reponce = 0, TR_min, TR_max, TR_top, u;
				u = Math.random();

				try {
					DEC = new TriangularFuzzyNumber(6, 10, 12);
					TR_min = -DEC.getA() * Math.log10(1 - u);
					TR_top = -DEC.getB() * Math.log10(1 - u);
					TR_max = -DEC.getC() * Math.log10(1 - u);
					MonteCarlo ran = new MonteCarlo();
					t_reponce = ran.simule(TR_min, TR_max, 2);
					System.out.println("T_responce= "+t_reponce);
				} catch (FuzzyOntologyException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				FirstOrderIntegrator dp853 = new DormandPrince853Integrator(
						6.14e-20, 1.0, 1.0e-10, 1.0e-10);
				FirstOrderDifferentialEquations ode = new NormalInputODE(Double
						.parseDouble(tbreak.getText()), Double.parseDouble(tmax
						.getText()), Double.parseDouble(k.getText()), Double
						.parseDouble(FcValue.getText()));
				
				double[] y = new double[] { Double.parseDouble(M.getText()),
						Double.parseDouble(Mc.getText()),
						Double.parseDouble(TR.getText()),
						Double.parseDouble(X.getText()),
						Double.parseDouble(open.getText()),
						Double.parseDouble(out.getText()) }; // initial state
				double[] ydot = new double[6];
				final double[][] tab = new double[1000][1];
				final double[][] tab1 = new double[1000][1];
				final double[][] tab2 = new double[1000][1];
				final double[][] tab3 = new double[1000][1];
				final double[][] tab4 = new double[1000][1];
				final double[][] tab5 = new double[1000][1];

				StepHandler stepHandler = new StepHandler() {

					ArrayList<String> steps = new ArrayList<String>();

					public void init(double t0, double[] y0, double t) {

					}

					public void handleStep(StepInterpolator interpolator,
							boolean isLast) {
						double t = interpolator.getCurrentTime();
						double[] y = interpolator.getInterpolatedState();

						steps.add(t + "-- " + y[0] + "-- " + y[1] + "-- "
								+ y[2] + "-- " + y[3] + " --" + y[4] + "-- "
								+ y[5]);

						somme += y[5];
				//		System.out.println(y[5]);
				//		System.out.println(y[2]);
						tab[(int) t][0] = y[0];
						tab1[(int) t][0] = y[1];
						tab2[(int) t][0] = y[2];
						tab3[(int) t][0] = y[3];
						tab4[(int) t][0] = y[4];
						tab5[(int) t][0] = y[5];
						if (isLast) {
							try {
								PrintWriter writer = new PrintWriter(
										new File(
												"C:/Users/houssein/Desktop/results1.txt"),
										"UTF-8");
								for (String step : steps) {
									writer.println(step);
								}
								writer.close();

							} catch (Exception e) {
							}
							;
						}

					}
				};

				dp853.addStepHandler(stepHandler);

				dp853.integrate(ode, 0.0, y, 900, y);

				massReleased.setText("mass released = " + somme);

				JavaPlot p = new JavaPlot();
				JavaPlot p1 = new JavaPlot();
				JavaPlot p2 = new JavaPlot();
				JavaPlot p3 = new JavaPlot();
				JavaPlot p4 = new JavaPlot();
				JavaPlot p5 = new JavaPlot();

				PlotStyle myPlotStyle = new PlotStyle();
				myPlotStyle.setStyle(Style.LINES);

				AbstractPlot s = new DataSetPlot(tab);
				AbstractPlot s1 = new DataSetPlot(tab1);
				AbstractPlot s2 = new DataSetPlot(tab2);
				AbstractPlot s3 = new DataSetPlot(tab3);
				AbstractPlot s4 = new DataSetPlot(tab4);
				AbstractPlot s5 = new DataSetPlot(tab5);

				myPlotStyle.setLineWidth(2);

				s.setTitle("Total mass");
				s1.setTitle("Oxide mass");
				s2.setTitle("Temperature");
				s3.setTitle("mass of oxide reacted");
				s4.setTitle("status of the burst disk");
				s5.setTitle("Vapor discharge rate");

				s.setPlotStyle(myPlotStyle);
				s1.setPlotStyle(myPlotStyle);
				s2.setPlotStyle(myPlotStyle);
				s3.setPlotStyle(myPlotStyle);
				s4.setPlotStyle(myPlotStyle);
				s5.setPlotStyle(myPlotStyle);
				p.addPlot(s);
				p1.addPlot(s1);
				p2.addPlot(s2);
				p3.addPlot(s3);
				p4.addPlot(s4);
				p5.addPlot(s5);
				p.newGraph();

				p1.newGraph();
				p2.newGraph();
				p3.newGraph();
				p4.newGraph();
				p5.newGraph();

				p.setTitle("Total mass");
				p1.setTitle("Oxide mass");
				p2.setTitle("Temperature");
				p3.setTitle("mass of oxide reacted");
				p4.setTitle("status of the burst disk");
				p5.setTitle("Vapor discharge rate");

				p.plot();
				p1.plot();
				p2.plot();
				p3.plot();
				p4.plot();
				p5.plot();

				M.setText("");
				Mc.setText("");
				TR.setText("");
				X.setText("");
				open.setText("");
				out.setText("");
				tbreak.setText("");
				tmax.setText("");
				FcValue.setText("");
				k.setText("");
			}
		});
		gl_contentPane.setHonorsVisibility(false);
		contentPane.setLayout(gl_contentPane);
	}
}
