package classes;
import org.apache.commons.math3.analysis.function.Exp;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;

import java.math.*;
import java.io.*;

public class NormalInputODE implements FirstOrderDifferentialEquations, EventHandler {
	public static double[][] points = new double[1000][1000];
	static int j = 1;
	double tbreak;
	double tmax;
	double FC;
	double K;

	public NormalInputODE(double tbreak, double tmax, double K, double FC) {
		this.tbreak = tbreak;
		this.tmax = tmax;
		this.K = K;
		this.FC = FC;
		

	}

	public int getDimension() {
		return 6;
	}

	public void computeDerivatives(double t, double[] y, double[] yDot) {
		// yDot[0] = omega * (c[1] - y[1]);
		// yDot[1] = omega * (y[0] - c[0]);
		double M = y[0], MC = y[1], TR = y[2], X = y[3], Open = y[4], out = y[5],

		Pburst = 9, T0 = 80, Lamda = 670, Cp = 3.5, P, V1, V, HR = -1660, R = 1.987, E = 21000, Kv = 100, M0 = 4400;

		double Fc = FC;
		if ((t > tbreak) & (t < tmax))
			Fc = 0;
	/*	if(t>705){
			Fc=5000;
		}*/
		// System.out.println(j+"."+Fc);
		double MW = (M0 + X) / (M0 / 74);

		double C = MC / M;

		double k = (9e9 * Math.exp(-E / (R * (TR + 273))))*K;
	
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
		double dTR = (Hc - Hv - Qg - Qr) / (M * Cp);
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

		// for (int i = 0; i < 6; i++) {
		// points[i][j] = yDot[i];
		// System.out.println(yDot[i]);
		// }
		j++;
		// System.out.println("************");
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
