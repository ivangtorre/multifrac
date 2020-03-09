// This method prints results after calculating the N(epsilon) in box counting monofractal method.

// Import libraries
import ij.plugin.filter.PlugInFilter;
import java.awt.*;
import ij.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.text.*;
import ij.gui.*;
import ij.util.*;
import ij.io.*;
import ij.process.*;
import ij.measure.*;
import static java.lang.Math.*;	

public class resultmonofractal{
double[] epsilon;
double[] N;
double w;
double h;
double maxsize;
boolean isstack;
boolean isgliding;
int numpixel;
String str;

	public void out() {
		// RESULTS TABLE
		// muestra: n, epsilon, (1/epsilon)^2, N(epsilon)
		// RESULTS TABLE
 		ResultsTable table = new ResultsTable();
 		for(int counter = 0; counter < epsilon.length; counter++)    {
			table.incrementCounter();
 			table.addValue("epsilon", epsilon[counter]);
			if (isstack == true && isgliding == false) {table.addValue("Total boxes", pow((1/epsilon[counter]),3));}
			else if (isstack == false && isgliding == false) {table.addValue("Total boxes", pow((1/epsilon[counter]),2));}
			else if (isstack == false && isgliding == true) {
			table.addValue("Size box", counter+1);
			table.addValue("Total boxes", (w-counter)*(h-counter));}
 			table.addValue("Occupied boxes", N[counter]);
		}
		table.show("Results for "+str);
		

		// CURVE FITTING (type of fitType: STRAIGHT_LINE=0, POLY2=1, POLY3=2, POLY4=3, EXPONENTIAL=4, POWER=5,...)	
		double[] x = new double[epsilon.length];					// x = log(1/epsilon)
		double[] y = new double[N.length];						// y = log(N)
		for(int counter = 0; counter < epsilon.length; counter++){
                	x[counter] = log(1/epsilon[counter]);
			y[counter] = log(N[counter]);
                }

		linefit fit1 = new linefit();
		double slope, devslope;
		CurveFitter cf;
		fit1.x   = x;									// slope parameter x
		fit1.y   = y;									// slope parameter y
		slope    = fit1.getslope();							// slope
		devslope = fit1.getdevslope();							// slope error
		cf       = fit1.getcf();							// contructor	
		plot(cf);									// plot

		// Print slope and its error
		IJ.log("");
		IJ.log("Box Counting Fractal Dimension:");
		IJ.log("");
		IJ.log("D = "+IJ.d2s(slope,5)+" +/- "+IJ.d2s(devslope,5));

	}

	public static void plot(CurveFitter cf) {
		plot(cf, false);
	}	

	// VARIABLES FOR THE PLOT	
	public static void plot(CurveFitter cf, boolean eightBitCalibrationPlot) {
		double[] x = cf.getXPoints();
		double[] y = cf.getYPoints();
		if (cf.getParams().length<cf.getNumParams()) {
			Plot plot = new Plot(cf.getFormula(),"ln(1/ε)","ln(N(ε))",x,y);
			plot.setColor(Color.BLUE);
			plot.addLabel(0.02, 0.1, cf.getName());
			plot.addLabel(0.02, 0.2, cf.getStatusString());
			plot.show();
			return;
		}
		// Trend line
		int npoints = 100;
		double[] a = Tools.getMinMax(x);
		double xmin=a[0], xmax=a[1];
		if (eightBitCalibrationPlot) {
			npoints = 256;
			xmin = 0;
			xmax = 255;
		}
		a = Tools.getMinMax(y);
		double ymin=a[0], ymax=a[1];
		float[] px = new float[npoints];
		float[] py = new float[npoints];
		double inc = (xmax-xmin)/(npoints-1);
		double tmp = xmin;
		for (int i=0; i<npoints; i++) {
			px[i]=(float)tmp;
			tmp += inc;
		}
		
		double[] params = cf.getParams();
		for (int i=0; i<npoints; i++)
			py[i] = (float)cf.f(params, px[i]);
		a    = Tools.getMinMax(py);
		ymin = min(ymin, a[0]);
		ymax = max(ymax, a[1]);
		Plot plot = new Plot("Fractal dimension","ln(1/ε)","ln(N(ε))",px,py);
		plot.setLimits(xmin, xmax, ymin, ymax);
		plot.setColor(Color.RED);
		plot.addPoints(x, y, PlotWindow.CIRCLE);
		plot.setColor(Color.BLUE);

		// Legend
		StringBuffer legend = new StringBuffer(100);
		legend.append(cf.getName()); legend.append('\n');
		legend.append(cf.getFormula()); legend.append('\n');
        	double[] p = cf.getParams();
        	int n = cf.getNumParams();
        	char pChar = 'a';
        	for (int i = 0; i < n; i++) {
			legend.append(pChar+" = "+IJ.d2s(p[i],5,9)+'\n');
			pChar++;
        	}
		legend.append("R^2 = "+IJ.d2s(cf.getRSquared(),4)); legend.append('\n');
		plot.addLabel(0.02, 0.1, legend.toString());
		plot.setColor(Color.BLUE);
		plot.show();									
	}



}

