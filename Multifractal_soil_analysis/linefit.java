// This class calculate the slope and its error (devslope) of two one-dimensional arrays
// Example of how it works:
// double[] x = {1,2,3,4,5};
// double[] y = {2,4,6,8,11};
// linefit fit1 = new linefit();
// fit1.x = x;
// fit1.y = y;
// double slope, devslope;
// slope = fit1.getslope();
// devslope = fit1.getdevslope();

// Import libraries
import ij.measure.*;
import static java.lang.Math.*;	

public class linefit {
	double[] x;
	double[] y;
	private double slope;
	private double origin;
	private double devslope;
	private static int fitType = 0;

	public CurveFitter getcf() {
		CurveFitter cf = new CurveFitter(x,y);
		cf.doFit(fitType);
		return cf;
	}

	public double getslope() {
		CurveFitter cf = new CurveFitter(x,y);
		cf.doFit(fitType);
		double[] aux = cf.getParams(); 							
		slope = aux[1];
		return slope;
	}

	public double getdevslope(){
		CurveFitter cf = new CurveFitter(x,y);
		cf.doFit(fitType);
		double[] aux = cf.getParams(); 							
		slope = aux[1];
		origin = aux[0];
		double Sx = 0;
		double Sy = 0;
		double Sxy = 0;
		double Sxx = 0;
		double chi=0;
		for (int i = 0; i <x.length; i++){
			Sx +=x[i];
			Sy +=y[i];
			Sxy +=x[i]*y[i];
			Sxx +=pow(x[i],2);
			chi += pow((y[i] - slope*x[i] - origin),2);
		}
		devslope = sqrt((x.length*chi)/((x.length*Sxx - Sx*Sx)*(x.length -2)));



		//double sigma;
		//sigma = cf.getSD();
		//double den1 = 0;
		//double den2 = 0;
		//for(int n = 0; n < x.length; n++){
		//	den1 += pow(x[n],2);
		//	den2 += x[n];
		//}
		//den1    *= x.length;
		//den2     = pow(den2,2);
		//devslope = (sqrt(x.length+1)*sigma)/sqrt(den1-den2);
		return devslope;
	}
	


}
				
