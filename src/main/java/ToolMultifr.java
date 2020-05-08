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

import ij.*;
import ij.gui.*;
import ij.measure.*;
import static java.lang.Math.*;

/**
 * Handles results for multifractal analysis
 * 
 * @author Ivan G Torre
 */

public class ToolMultifr{
	double[] epsilon;
	double[][] Xnum;
	double[][] denomMu;  
	double[][] numalpha;  
	double[][] numf;
	double[] q;
	boolean isstack;
	String str;
	String path;

	public void out() {
		double[] x = new double[epsilon.length];
		
		double minx0 = 0;
		double maxx0 = 0;
		for(int n = 0; n < epsilon.length; n++){
			x[n] = log(epsilon[n]);
			if (x[n]>maxx0){maxx0=x[n];}
			if (x[n]<minx0){minx0=x[n];}
		}
		
		double miny0 = 0;
		double maxy0 = 0;
		double[][] y = new double[epsilon.length][q.length];
		for (int n = 0; n<epsilon.length; n++){
			for (int qposition = 0; qposition<q.length; qposition++){
				if (q[qposition] != 1){
					y[n][qposition] = log(Xnum[n][qposition]);
				}
				else{	
					y[n][qposition] += Xnum[n][qposition];
				}
				if (y[n][qposition]>maxy0){maxy0=y[n][qposition];}
				else if (y[n][qposition]<miny0){miny0=y[n][qposition];}	 	
			}
		}			

		
    	/** ---------Compute Multifractal Features---------------------------------- */         

		double[] TauQ     	= new double[q.length];
		double[] desvTauQ 	= new double[q.length];
		double[] Dq       	= new double[q.length];
		double[] desvDq   	= new double[q.length];
		double[] alphaq   	= new double[q.length];
		double[] desvalphaq   	= new double[q.length];
		double[] fq  		= new double[q.length];
		double[] desvfq   	= new double[q.length];
		double[] ybar = new double[y.length];

		/** ---------Tauq---------------------------------- */ 		
		for(int qposition = 0; qposition < q.length; qposition++){
			double[] yaux = new double[y.length];
			for(int n = 0; n < y.length; n++){
				yaux[n] = y[n][qposition];
			}
			
			Linefit fit1 = new Linefit();
			fit1.x   = x;
			fit1.y   = yaux;
			TauQ[qposition]     = fit1.getslope();
			desvTauQ[qposition] = fit1.getdevslope();	  
		}

		/** ---------Dq---------------------------------- */ 		
		for (int qposition = 0; qposition<q.length; qposition++){
			if (q[qposition] != 1){
				Dq[qposition]     = TauQ[qposition]/(q[qposition] - 1);
				desvDq[qposition] = desvTauQ[qposition]/abs(q[qposition] - 1);
			}
			else if (q[qposition] == 1){
				Dq[qposition]       = TauQ[qposition];
				desvDq[qposition]   = desvTauQ[qposition];
				TauQ[qposition]     = 0;
			}	
		}
		
		/** ---------alfa---------------------------------- */ 		
		double maxx = 0;
		double minx = 10;
		double maxdesvalphaq = 0;
	
		for(int qposition = 0; qposition < q.length; qposition++){
			double[] yalphaq = new double[epsilon.length];
			for(int n = 0; n < epsilon.length; n++){
				yalphaq[n] = numalpha[n][qposition];
			}
			Linefit fit2 = new Linefit();
			fit2.x   = x;
			fit2.y   = yalphaq;
			alphaq[qposition]    = fit2.getslope();
			if (alphaq[qposition]<minx){minx=alphaq[qposition];}
			if (alphaq[qposition]>maxx){maxx=alphaq[qposition];}
			desvalphaq[qposition] = fit2.getdevslope();
		  	if (desvalphaq[qposition]>maxdesvalphaq){maxdesvalphaq=desvalphaq[qposition];}
		}

		/** ---------f(alfa)---------------------------------- */ 		
		double maxy = 0;
		double miny = 10;
		double maxdesvfq = 0;
		for(int qposition = 0; qposition < q.length; qposition++){
			double[] yfq = new double[epsilon.length];
			for(int n = 0; n < epsilon.length; n++){
				yfq[n] = numf[n][qposition];
			}
			Linefit fit3 = new Linefit();
			fit3.x   = x;
			fit3.y   = yfq;
			fq[qposition]    = fit3.getslope();
			if (fq[qposition]<miny){miny=fq[qposition];}
			if (fq[qposition]>maxy){maxy=fq[qposition];}
			desvfq[qposition] = fit3.getdevslope();
			if (desvfq[qposition]>maxdesvfq){maxdesvfq=desvfq[qposition];}
		}

		
		/** --------PLOTS---------------------------------- */ 		
		/** --------Plot Xq-epsilon---------------------------------- */ 		
		Plot plot0 = new Plot("X(q,epsilon)"+str, "ln(epsilon)", "ln(X(q))");
		PlotWindow.noGridLines = false;
		plot0.setLimits(minx0-0.5, maxx0 + 0.5,miny0-1,maxy0+1);
		for(int qposition = 0; qposition < q.length; qposition++){
			for(int n = 0; n < y.length; n++){
				ybar[n] = y[n][qposition];
			}
		plot0.addPoints(x, ybar, Plot.LINE);
		}
		plot0.show();
		ImagePlus imp_plot0 = plot0.getImagePlus();
		IJ.saveAs(imp_plot0, "tif", path + "X(q,epsilon)"+ str);

		/** --------Plot Tauq---------------------------------- */ 		
		PlotWindow.noGridLines = false;
		Plot plot1 = new Plot("Tau(q)"+str, "q", "TauQ");
		plot1.add("line", q, TauQ);
		plot1.setLimits((q[0])-0.5, q[q.length -1]+0.5, TauQ[0]-desvTauQ[0]-0.5, TauQ[q.length -1]+desvTauQ[0]+0.5);
		plot1.setLineWidth(2);
		plot1.addErrorBars(desvTauQ); 
		plot1.show();
		ImagePlus imp_plot1 = plot1.getImagePlus();
		IJ.saveAs(imp_plot1, "tif", path + "TauQ vs q"+ str);

		/** --------Plot Dq---------------------------------- */ 		
		PlotWindow.noGridLines = false;
		Plot plot2 = new Plot("D(q)"+str, "q", "Dq");
		plot2.add("line",q, Dq);
		plot2.setLimits((q[0]) -0.5, q[q.length -1] + 0.5, Dq[q.length -1] - desvDq[q.length -1] - 0.5, Dq[0] + desvDq[0] + 0.5);
		plot2.setLineWidth(2);
		plot2.addErrorBars(desvDq);
		plot2.show();
		ImagePlus imp_plot2 = plot2.getImagePlus();
		IJ.saveAs(imp_plot2, "tif", path + "Dq vs q"+ str);

		/** --------Plot f(alfa)---------------------------------- */ 		
		PlotWindow.noGridLines = false;
		Plot plot3 = new Plot("f(alpha) vs alpha"+str, "alpha", "f(alpha)");
		plot3.add("line", alphaq, fq);
		plot3.setLimits(minx-maxdesvalphaq,maxx+maxdesvalphaq,miny-maxdesvfq,maxy+maxdesvfq);
		plot3.setLineWidth(2);
		plot3.addErrorBars(desvalphaq);
		for (int qposition = 0; qposition<q.length; qposition++){
			plot3.drawLine(alphaq[qposition]-desvalphaq[qposition],fq[qposition],alphaq[qposition]+desvalphaq[qposition],fq[qposition]);
		}
		plot3.show();
		ImagePlus imp_plot3 = plot3.getImagePlus();
		IJ.saveAs(imp_plot3, "tif", path + "f(alpha) vs alpha"+ str);

		/** --------Results Table---------------------------------- */ 		
 		ResultsTable table1 = new ResultsTable();
 		for(int counter = 0; counter < q.length; counter++)    {
			table1.incrementCounter();
			table1.addValue("q", IJ.d2s(q[counter], 2));
			table1.addValue("TauQ",IJ.d2s(TauQ[counter], 8));
			table1.addValue("error TauQ",IJ.d2s(desvTauQ[counter],8));
			table1.addValue("Dq", IJ.d2s(Dq[counter],8));
			table1.addValue("error Dq", IJ.d2s(desvDq[counter], 8));
 			table1.addValue("alpha", IJ.d2s(alphaq[counter], 8));
 			table1.addValue("error alpha", IJ.d2s(desvalphaq[counter], 8));
 			table1.addValue("f(alpha)", IJ.d2s(fq[counter], 8));
			table1.addValue("error f", IJ.d2s(desvfq[counter],8));
		}

		for(int ep = 0; ep < epsilon.length; ep++){
			table1.incrementCounter();
			for(int n = 0; n < q.length; n++){
				String fs;
				fs = String.format("x %d",n);
				table1.addValue(fs, x[ep]);
				String fs2;
				fs2 = String.format("y %d",n);
				table1.addValue(fs2, y[ep][n]);
					
			}
		}
		table1.setPrecision(8);
		table1.updateResults();
		table1.show("Results");
		IJ.selectWindow("Results");
		IJ.saveAs("Results", path + "Results.csv");
	}
}
