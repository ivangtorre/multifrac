/*
#%L
MULTIFRAC
ImageJ Plugin for fractal, multifractal and scaling image characterization and analysis.
Copyright (C) 2020  Ivan G Torre
-------Cite:-----------------------------------------------------------------------
I. G. Torre and A. M. Tarquis, MULTIFRAC: An ImageJ plugin for multifractal and 
multiscale characterization of 2D and 3D stack images, preprint.
------------------------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/
import java.awt.*;
import ij.*;
import ij.gui.*;
import ij.util.*;
import ij.measure.*;
import static java.lang.Math.*;

/**
 * Tool that handle results from Fractal Dimensions
 * 
 * @author Ivan G Torre
 */

public class ToolFracDim{
double[] epsilon, N;
double w, h,  maxsize;
boolean isstack, isgliding;
int numpixel;
static String str;
static String path;

	public void out() {
		
		/** ----------Table Results ----------------------------- */
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
		table.setPrecision(8);
		table.updateResults();
		table.show("Results");
		IJ.selectWindow("Results");
		IJ.saveAs("Results", path + "Results.csv");

		
		/** ----------Fit results line----------------------------- */
		double[] x = new double[epsilon.length];
		double[] y = new double[N.length];
		for(int counter = 0; counter < epsilon.length; counter++){
                	x[counter] = log(1/epsilon[counter]);
			y[counter] = log(N[counter]);
                }

		Linefit fit1 = new Linefit();
		double slope, devslope;
		CurveFitter cf;
		fit1.x   = x;
		fit1.y   = y;
		slope    = fit1.getslope();
		devslope = fit1.getdevslope();
		cf       = fit1.getcf();
		plot(cf);
		IJ.log("Box Counting Fractal Dimension:");
		IJ.log("D = "+IJ.d2s(slope,5)+" +/- "+IJ.d2s(devslope,5));

	}

	public static void plot(CurveFitter cf) {
		plot(cf, false);
	}	

	/** ----------Variables for the plot--------------------------------- */
	public static void plot(CurveFitter cf, boolean eightBitCalibrationPlot) {
		double[] x = cf.getXPoints();
		double[] y = cf.getYPoints();
		if (cf.getParams().length<cf.getNumParams()) {
			Plot plot = new Plot(cf.getFormula(), "ln(1/ε)", "ln(N(ε))");
			plot.add("line", x, y);	
			plot.addLabel(0.02, 0.1, cf.getName());
			plot.addLabel(0.02, 0.2, cf.getStatusString());
			plot.show();
			return;
		}
		/** ----------Trend--------------------------------- */
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
		double[] px = new double[npoints];
		double[] py = new double[npoints];
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
		Plot plot = new Plot("Fractal dimension", "ln(1/ε)", "ln(N(ε))");
		plot.add("line",px, py);
		plot.setLimits(xmin, xmax, ymin, ymax);
		plot.setColor(Color.RED);
		plot.add("box", x, y);
		plot.setColor(Color.BLUE);

		/** ----------Legend--------------------------------- */
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
		
		/** ----------Save plot--------------------------------- */
		ImagePlus imp_plot = plot.getImagePlus();
		IJ.saveAs(imp_plot, "tif", path + "Fractal_dimension"+ str);
	}

}

