// This method prints results after calculating the N(epsilon) in box counting monofractal method.

// Import libraries
import ij.plugin.filter.PlugInFilter;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;
import ij.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.text.*;
import ij.gui.*;
import ij.util.*;
import ij.io.*;
import ij.process.*;
import ij.measure.*;
import javax.imageio.ImageIO;
import static java.lang.Math.*;	

public class resultmultifractal{
	double[] epsilon;
	double[][] Xnum;
	double[][] denomMu;  
	double[][] numalpha;  
	double[][] numf;
	double[] q;
	boolean isstack;
	String str;
	String path;

	// Final variables x e y
	public void out() {
		double[] x = new double[epsilon.length];					// x = log(epsilon)
		
		double minx0 = 0;
		double maxx0 = 0;
		for(int n = 0; n < epsilon.length; n++){
			x[n] = log(epsilon[n]);
			if (x[n]>maxx0){maxx0=x[n];}
			if (x[n]<minx0){minx0=x[n];}
		}
		
		double miny0 = 0;
		double maxy0 = 0;
		double[][] y = new double[epsilon.length][q.length]; 				// y = log(X(q))
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


		// CURVE FITTING (type of fitType: STRAIGHT_LINE=0, POLY2=1, POLY3=2, POLY4=3, EXPONENTIAL=4, POWER=5,...)	
		double[] TauQ     	= new double[q.length];					// Slopes
		double[] desvTauQ 	= new double[q.length];					// Tau error
		double[] Dq       	= new double[q.length];					// Dq
		double[] desvDq   	= new double[q.length];					// Dq error
		double[] alphaq   	= new double[q.length];					// Dq error
		double[] desvalphaq   	= new double[q.length];					// Dq error
		double[] fq  		= new double[q.length];					// Dq error
		double[] desvfq   	= new double[q.length];					// Dq error
			
		double[] ybar = new double[y.length];

		
		// CALCULATIONS //

		// TauQ and eror
		for(int qposition = 0; qposition < q.length; qposition++){
			double[] yaux = new double[y.length];					// Adjusting aux array
			for(int n = 0; n < y.length; n++){
				yaux[n] = y[n][qposition];
			}
			

			linefit fit1 = new linefit();
			fit1.x   = x;								// slope parameter x
			fit1.y   = yaux;							// slope parameter y
			TauQ[qposition]     = fit1.getslope();					// slope
			desvTauQ[qposition] = fit1.getdevslope();				// slope error		  
		}

		// Dq and error
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
		
		// Alhpa (q) and error
		double maxx = 0;
		double minx = 10;
		double maxdesvalphaq = 0;
	
		for(int qposition = 0; qposition < q.length; qposition++){
			double[] yalphaq = new double[epsilon.length];				// Adjusting aux array
			for(int n = 0; n < epsilon.length; n++){
				yalphaq[n] = numalpha[n][qposition];
			}
			linefit fit2 = new linefit();
			fit2.x   = x;								// slope parameter x
			fit2.y   = yalphaq;							// slope parameter y
			alphaq[qposition]    = fit2.getslope();					// slope
			if (alphaq[qposition]<minx){minx=alphaq[qposition];}
			if (alphaq[qposition]>maxx){maxx=alphaq[qposition];}
			desvalphaq[qposition] = fit2.getdevslope();				//slope error
		  	if (desvalphaq[qposition]>maxdesvalphaq){maxdesvalphaq=desvalphaq[qposition];}
		}

		// f (q) and error
		double maxy = 0;
		double miny = 10;
		double maxdesvfq = 0;
		for(int qposition = 0; qposition < q.length; qposition++){
			double[] yfq = new double[epsilon.length];				// Adjusting aux array
			for(int n = 0; n < epsilon.length; n++){
				yfq[n] = numf[n][qposition];
			}
			linefit fit3 = new linefit();
			fit3.x   = x;								// slope parameter x
			fit3.y   = yfq;								// slope parameter y
			fq[qposition]    = fit3.getslope();					// slope
			if (fq[qposition]<miny){miny=fq[qposition];}
			if (fq[qposition]>maxy){maxy=fq[qposition];}
			desvfq[qposition] = fit3.getdevslope();					// slope error  
			if (desvfq[qposition]>maxdesvfq){maxdesvfq=desvfq[qposition];}
		}


		// PLOTS //

		// Plot Xq vs epsilon
		Plot plot0 = new Plot("X(q,epsilon)"+str, "ln(epsilon)", "ln(X(q))");
		PlotWindow.noGridLines = false; 		// draw grid lines
		plot0.setLimits(minx0-0.5, maxx0 + 0.5,miny0-1,maxy0+1);
		for(int qposition = 0; qposition < q.length; qposition++){
			for(int n = 0; n < y.length; n++){
				ybar[n] = y[n][qposition];
			}
		plot0.addPoints(x, ybar, Plot.LINE);
		}
		plot0.show();
		ImagePlus imp_plot0 = plot0.getImagePlus();
		IJ.saveAs(imp_plot0, "Jpeg", path + "X(q,epsilon)"+ str);
		imp_plot0.close();
		IJ.selectWindow("X(q,epsilon)"+str);




		// Plot TauQ vs q
		PlotWindow.noGridLines = false; 						// draw grid lines
		Plot plot1 = new Plot("TauQ vs q"+str, "q", "TauQ", q, TauQ);
		plot1.setLimits((q[0])-0.5, q[q.length -1]+0.5, TauQ[0]-desvTauQ[0]-0.5, TauQ[q.length -1]+desvTauQ[0]+0.5);
		plot1.setLineWidth(2);								// Line width
		plot1.addErrorBars(desvTauQ); 							// Errors bars
		plot1.setColor(Color.red);							// line color
		plot1.show();
		ImagePlus imp_plot1 = plot1.getImagePlus();
		IJ.saveAs(imp_plot1, "Jpeg", path + "TauQ vs q"+ str);
		imp_plot1.close();		



		// Plot q vs Dq
		PlotWindow.noGridLines = false; 						// draw grid lines
		Plot plot2 = new Plot("Dq vs q"+str, "q", "Dq", q, Dq);
		plot2.setLimits((q[0]) -0.5, q[q.length -1] + 0.5, Dq[q.length -1] - desvDq[q.length -1] - 0.5, Dq[0] + desvDq[0] + 0.5);
		plot2.setLineWidth(2);								// Line width
		plot2.addErrorBars(desvDq); 							// Errors bars
		plot2.setColor(Color.red);							// line color
		plot2.show();
		ImagePlus imp_plot2 = plot2.getImagePlus();
		IJ.saveAs(imp_plot2, "Jpeg", path + "Dq vs q"+ str);
		imp_plot2.close();


		// Plot f vs alpha
		PlotWindow.noGridLines = false; 						// draw grid lines
		Plot plot3 = new Plot("f(alpha) vs alpha"+str, "alpha", "f(alpha)", alphaq, fq);
		plot3.setLimits(minx-maxdesvalphaq,maxx+maxdesvalphaq,miny-maxdesvfq,maxy+maxdesvfq);
		plot3.setLineWidth(2);								// Line width
		plot3.addErrorBars(desvalphaq); 						// Errors bars
		plot3.setColor(Color.red);							// line color
		for (int qposition = 0; qposition<q.length; qposition++){
			plot3.drawLine(alphaq[qposition]-desvalphaq[qposition],fq[qposition],alphaq[qposition]+desvalphaq[qposition],fq[qposition]);
		}
		plot3.show();
		ImagePlus imp_plot3 = plot3.getImagePlus();
		IJ.saveAs(imp_plot3, "Jpeg", path + "f(alpha) vs alpha"+ str);
		imp_plot3.close();

		// RESULTS TABLE
		// muestra: alpha, error alpha, f(alpha), error f(alpha)
		// RESULTS TABLE
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
		table1.setPrecision(8);


		table1.show("Results");
		IJ.selectWindow("Results");
		IJ.save(path+"Results.csv");
		




	}
}
