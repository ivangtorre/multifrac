/*
This is a method from _2D_ Structure_Function and _3D_Structure Function
programs. This calculate the results and plot the struture functions.
The requiriments are Mq, lista_tau and lista_q.
Programed by Ivan Gonzalez Torre, ivangonzaleztorre@gmail.com
*/

// Libraries import
import java.awt.*;
import ij.gui.*;
import static java.lang.Math.*;

public class resultestructure{

	double[][] Mq;
	int[] lista_tau;
	double[] lista_q;

	public void out() {
	
		// ZETA calculation for each q
		double[] Zeta 		= new double [lista_q.length];
		double[] Hurst 		= new double [lista_q.length];
		double[] desvZeta 	= new double [lista_q.length];
		double[] desvHurst 	= new double [lista_q.length];		
		double[] xaux 		= new double [lista_tau.length];
		double minx = 0, miny = 20, maxy = -20, minZeta = 0, maxZeta = 0, maxdesvZeta = 0, minHurst = 10, maxHurst = 0;

		// Values of X
		for(int i = 0; i < lista_tau.length; i++){
			xaux[i] = log((double)lista_tau[i]/lista_tau[lista_tau.length - 1]);
			if (xaux[i]<minx){minx=xaux[i];}			// Axis limits for plotting
		}
		
		// Values of Y
		for (int qposition = 0; qposition < lista_q.length; qposition++){
			double[] yaux = new double [lista_tau.length];
			for(int i = 0; i < lista_tau.length; i++){
				yaux[i] = log(Mq[i][qposition]);
				if (yaux[i]>maxy){maxy=yaux[i];}		// Axis limits for plotting
				if (yaux[i]<miny){miny=yaux[i];}		// Axis limits for plotting
			}
			// Fit the slope
			linefit fit1 = new linefit();
			fit1.x   = xaux;					
			fit1.y   = yaux;					
			Zeta[qposition]     = fit1.getslope();					// slope
			if (Zeta[qposition]>maxZeta){maxZeta=Zeta[qposition];}
			if (Zeta[qposition]/lista_q[qposition]>maxHurst){maxHurst=Zeta[qposition]/lista_q[qposition];}
			if (Zeta[qposition]<minZeta){minZeta=Zeta[qposition];}
			if (Zeta[qposition]/lista_q[qposition]<minHurst){minHurst=Zeta[qposition]/lista_q[qposition];}
			desvZeta[qposition] = fit1.getdevslope();				// slope error
			if (desvZeta[qposition]>maxdesvZeta){maxdesvZeta=desvZeta[qposition];}
		}
		

		// Hurst exponent
		// The errors of Hurst exponent are the same that the Zeta function
		for (int qposition = 0; qposition < lista_q.length; qposition++){
			Hurst[qposition] = Zeta[qposition]/lista_q[qposition];
			desvHurst[qposition] = desvZeta[qposition]/lista_q[qposition];
		}


		////////////////////////////////////////////////////
		///////////////////// PLOT /////////////////////////
		////////////////////////////////////////////////////
		// Plot log(Tau/Taumax) vs log(Mq(Tau))
		Plot plot0 = new Plot("Mq(tau)", "ln(Tau/Taux_max)", "ln(Mq(Tau))");
		PlotWindow.noGridLines = false; 					// draw grid lines
		plot0.setLimits(minx-0.5, 0, miny-0.5, maxy+0.5);
		for (int qposition = 0; qposition < lista_q.length; qposition++){
			double[] yaux = new double [lista_tau.length];
			for(int i = 0; i < lista_tau.length; i++){
				yaux[i] = log(Mq[i][qposition]);
			}
			plot0.addPoints(xaux, yaux, Plot.LINE);
		}

		plot0.show();

  		// PLOT ZETA VS Q
		// Line y=0.5x
		double[] y05 = new double [2];
		double[] x05 = new double [2];
		x05[0] = 0;
		x05[1] = lista_q[lista_q.length -1];
		y05[0] = 0;
		y05[1] = 0.5*x05[1];

		// Line from (0,0) to (f(q=1),q=1)
		double[] yq1 = new double [2];
		double[] xq1 = new double [2];
		xq1[0] = 0;
		xq1[1] = lista_q[lista_q.length -1];
		yq1[0] = 0;
		double slope_yq1=0;
		
		
		// This loop looks the value q=1
		boolean isq1 = false; // Check is exist q=1, if not interpolate
		int position = 0;
		for (int qposition = 0; qposition < lista_q.length; qposition++){
			if (lista_q[qposition] == 1){
				position = qposition;
				isq1 = true;
				slope_yq1 = Zeta[position];
			}
		}
		// If not exist the value q =1, interpolate it.
		if (isq1 == false){
			int pos1 = 0, pos2 = 0;
			double slope_aux;
			for (int qposition = 0; qposition < lista_q.length; qposition++){
				if (lista_q[qposition] <= 1){
					pos1 = qposition;
				}
			}
			pos2 = pos1 + 1;
			slope_aux = (Zeta[pos2] - Zeta[pos1])/(lista_q[pos2] - lista_q[pos1]);
			slope_yq1 = Zeta[pos1] + slope_aux*(1-lista_q[pos1]);
		}

		yq1[1] = slope_yq1*xq1[1];
		
		//Plot
		PlotWindow.noGridLines = false; 						// draw grid lines
		Plot plot1 = new Plot("Structure function", "q", "Zeta", lista_q, Zeta);
		plot1.setLimits((lista_q[0])-0.5, lista_q[lista_q.length -1]+0.5, minZeta-0.5-maxdesvZeta, maxZeta+0.5+maxdesvZeta);
		plot1.setLineWidth(1);								// Line width
		plot1.addErrorBars(desvZeta); 							// Errors bars
		plot1.setColor(Color.magenta);							// line color
		plot1.addPoints(x05, y05, Plot.LINE);
		plot1.setColor(Color.blue);
		plot1.addPoints(xq1, yq1, Plot.LINE);		
		plot1.setColor(Color.red);
		plot1.setLineWidth(2);
		double xloc = 0.7;
		double yloc = 0.8;
		plot1.setColor(Color.black);
		plot1.setJustification(Plot.RIGHT);
		plot1.addLabel(xloc, yloc, "Structure Function");
		plot1.addLabel(xloc, yloc + 0.06,"y = 0.5*x");
		plot1.addLabel(xloc, yloc + 0.12,"y = (Zeta(q=1))*x");
		xloc += 0.01;
		yloc -= 0.01;
		plot1.setColor(Color.red);
		plot1.drawNormalizedLine(xloc, yloc-0.02, xloc+0.1, yloc-0.02);
		plot1.setColor(Color.magenta);
		plot1.drawNormalizedLine(xloc, yloc-0.02+0.06, xloc+0.1, yloc-0.02+0.06);
		plot1.setColor(Color.blue);
		plot1.drawNormalizedLine(xloc, yloc-0.02+0.12, xloc+0.1, yloc-0.02+0.12);

		plot1.setColor(Color.red);
		plot1.show();


		// Plot Hurst
		PlotWindow.noGridLines = false; 						// draw grid lines
		Plot plot2 = new Plot("Hurst exponents", "q", "H(q)", lista_q, Hurst);
		plot2.setLimits((lista_q[0])-0.5, lista_q[lista_q.length -1]+0.5, minHurst-0.3-maxdesvZeta, maxHurst+0.3);
		plot2.setLineWidth(1);								// Line width
		plot2.addErrorBars(desvHurst); 							// Errors bars
		plot2.setLineWidth(2);
		plot2.show();	
	
		


	}
}
