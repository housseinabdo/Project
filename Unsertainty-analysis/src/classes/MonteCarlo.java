package classes;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class MonteCarlo {

/**
* @param args
*/
public Double simule(double p1, double p2, int loi) {
Random randomno = new Random();
double y = 0.0;

if (loi == 1) {
// loi de distrubtion uniforme

y = randomno.nextDouble() * (p2 - p1) + p1;
} else if (loi == 2) {
// loi normale

double Moy = (p1 + p2) / 2;
double var = ((p2 - p1) / 6) * ((p2 - p1) / 6);

y = randomno.nextGaussian() * Math.sqrt(var) + Moy;
}

return y;
}

}