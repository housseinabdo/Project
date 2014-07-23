package classes;
import org.apache.commons.math3.analysis.function.Exp;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;

import fuzzydl.TriangularFuzzyNumber;
import fuzzydl.exception.FuzzyOntologyException;

import java.math.*;
import java.io.*;

public class CircleODE implements FirstOrderDifferentialEquations, EventHandler {
	public static double[][] points = new double[1000][1000];
	double tbreak , t_reponce;
	double f_deluge;
	double tmax;
	double FC;
	double K;
	TriangularFuzzyNumber DEC;

	public CircleODE(double f_deluge, double K, double FC,double t_break,double t_reponce) {
		this.tbreak=t_break;
		this.t_reponce=t_reponce;
		this.f_deluge = f_deluge;
		this.K = K;
		this.FC = FC;
		
	}
	public CircleODE(double tbreak, double tmax, double K, double FC) {

		this.tbreak = tbreak;
		this.tmax = tmax;
		this.K = K;
		this.FC = FC;

	}
	public int getDimension() {
		return 6;
	}

	public void computeDerivatives(double t, double[] y, double[] yDot) {
		double M = y[0], MC = y[1], TR = y[2], X = y[3], Open = y[4], out = y[5], F_deluge, A = 1,

		Pburst = 9, T0 = 80, Lamda = 670, Cp = 3.5, P, V1, V, HR = -1660, R = 1.987, E = 21000, Kv = 100, M0 = 4400;

		double Fc = FC, Q_deluge = 0;
		if ((t > tbreak)) {
			Fc = 0;
		}
		if (t > tbreak + t_reponce) {
			F_deluge = this.f_deluge;
			Q_deluge = F_deluge * A * Cp * (TR - T0);
		}

		double MW = (M0 + X) / (M0 / 74);

		double C = MC / M;

		double k = (9e9 * Math.exp(-E / (R * (TR + 273)))) * K;

		double P1 = (Math.exp(-3430 / (TR + 273) + 11.7) + 1.45e-3 * MW) * C;

		if (P1 < 1)
			P = 1;
		else
			P = P1;

		double Qr = Fc * Cp * (TR - T0);

		double r = k * MC;
		double Qg = r * HR;

		double Vsubs = Kv * P / Math.sqrt((TR + 273))
				* Math.sqrt(1 + 1 / (P * P));
		double Vs = 0.85 * Kv * P / Math.sqrt(TR + 273);

		if (P < 1.9)
			V1 = Vsubs;
		else
			V1 = Vs;

		if ((P <= 1) | (Open == 0))
			V = 0;
		else
			V = V1;

		double F, dOpen;
		if (Open > 0)
			F = 0;

		else
			F = 100;

		double Hv = V * Lamda;
		double Hc = F * Cp * (T0 - TR);

		if (P < Pburst)
			dOpen = 0;
		else
			dOpen = 0.001;

		double dX = r;
		double dTR = (Hc - Hv - Qg - Qr - Q_deluge) / (M * Cp);
		double dMC = F - V - r;
		double dM = F - V;

		// Relation ajoutée
		double dOut = V;
		yDot[0] = dM;
		yDot[1] = dMC;
		yDot[2] = dTR;
		yDot[3] = dX;
		yDot[4] = dOpen;
		yDot[5] = V;

	}
	
	@Override
	public Action eventOccurred(double arg0, double[] arg1, boolean arg2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double g(double arg0, double[] arg1) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void init(double arg0, double[] arg1, double arg2) {
		// TODO Auto-generated method stub

	}

	@Override
	public void resetState(double arg0, double[] arg1) {
		// TODO Auto-generated method stub

	}

}
