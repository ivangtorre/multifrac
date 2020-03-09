// Ivan Gonzalez Torre

// Libraries import
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
import java.util.Arrays;


public class _3D_Gliding_Multifractal implements PlugInFilter {

protected ImageStack stackaux;

	public int setup(String arg, ImagePlus impaux) {
		stackaux = impaux.getStack();
		if (impaux.getProcessor().isInvertedLut()){
			IJ.run("Invert LUT");
		}
		return DOES_8G+DOES_16+DOES_32+STACK_REQUIRED;				// Suports gray scale image
	}


	public void run(ImageProcessor ip) {						// Process image and save in IP
		
		String path = "/home/ivangt/Desktop/";

		int resolution = ip.getBitDepth();
		// Input dialog for q values
		GenericDialog gdq = new GenericDialog("Select q values for multifractal dimensions");	// Create dialog
		gdq.addNumericField("qmin: ", -10, 0);							// Field1
		gdq.addNumericField("qmax: ", 10, 0);							// Field2
		gdq.addNumericField("qincrement: ", 0.1, 1);						// Field3
		gdq.showDialog();									// Show
		if (gdq.wasCanceled()) return;								// If no input, cancell program
		double qmin, qmax, qincrement;
		qmin       = (double)gdq.getNextNumber();						// Save value1
		qmax       = (double)gdq.getNextNumber();						// Save value2
		qincrement = (double)gdq.getNextNumber();						// Save value3
		double[] q = new double[(int)((qmax - qmin)/qincrement)+1];				// Declare q array
		for (int count = 0; count<q.length; count++){						// Asign q values to the array	
			q[count] = qmin + count*qincrement;
		}
			
		// Check size of image
		double woriginal = ip.getWidth();					// Width from the original image
		double horiginal  = ip.getHeight();					// heigth from the original image
		double doriginal = stackaux.getSize();					// depth from the original image


		// Creating the new image
		IJ.run("Options...", "iterations=1 black count=1"); 
		IJ.run("Duplicate...","title=Analyzed_image.tif duplicate");	// Duplicate image
		IJ.run("Options...", "iterations=1 black count=1");			// Check not invertign LUT

		
		GenericDialog gd1	= new GenericDialog("Process image");			// Create dialog
		String[] titleArray2 	= new String[2];
		gd1.addMessage("Maximum measure will be for the black. Do you want to analyze this image or the inverted one?");
		titleArray2[0] = "This image";
		titleArray2[1] = "Inverted image";
		gd1.addChoice("  ", titleArray2, titleArray2[0]);
		gd1.showDialog();
		if (gd1.wasCanceled()) return;								// If no input, cancell program
		int choose;
		choose = gd1.getNextChoiceIndex();
		if (choose == 0){
		}
		else if (choose == 1){// Then invert image
			ImagePlus imp = IJ.getImage(); 
			ImageStack stack = imp.getStack(); 
			Roi cuadradointeres;							// variable tipo Roi
			Rectangle roi;								// Region of interest multifractal
			ImageProcessor salida_ip = stack.getProcessor(1);
			Rectangle fullimage = salida_ip.getRoi();

			if (resolution == 16){
				Rectangle r = salida_ip.getRoi();
				for (int z=1; z<stack.getSize()+1; z++){
					salida_ip = stack.getProcessor(z);
	      				for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){// Normalize matriz
						for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
		  					salida_ip.set(x, y, (int)abs(65535 - salida_ip.get(x,y)));
	    					}
	  				}
				}
			imp.updateAndDraw();
			IJ.wait(500);			
			} 

			else if (resolution != 16){
				IJ.run("Invert", "stack");				// Invert image. Count black mass
				IJ.wait(500);
			}
		}

		ImagePlus imp = IJ.getImage(); 
		ImageStack stack = imp.getStack(); 
		Roi cuadradointeres;							// variable tipo Roi
		Rectangle roi;								// Region of interest multifractal
		ImageProcessor salida_ip = stack.getProcessor(1);
		Rectangle fullimage = salida_ip.getRoi();

		if (resolution == 16){
			Rectangle r = salida_ip.getRoi();
			for (int z=1; z<stack.getSize()+1; z++){
				salida_ip = stack.getProcessor(z);
      				for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){// Normalize matriz
					for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
          					salida_ip.set(x, y, (int)abs(65535 - salida_ip.get(x,y)));
    					}
  				}
			}
		imp.updateAndDraw();
		IJ.wait(500);
		} 

		else if (resolution != 16){
			IJ.run("Invert", "stack");				// Invert image. Count black mass
			IJ.wait(500);
		}

		long time_start, time_end;
		time_start = System.currentTimeMillis();

		// Begining of the program
		double sum = 0;								// Sum of all pixels value 
		for (int z=1; z<stack.getSize()+1; z++){ 
			ImageProcessor salidaux_ip = stack.getProcessor(z);
			sum += salidaux_ip.getStatistics().area*salidaux_ip.getStatistics().mean;
		}
		
		double maxsize	    	= min(woriginal,horiginal);				
		maxsize	    		= min(woriginal,doriginal);			// Max size of ROI
		double[] epsilon 	= new double[(int)maxsize+1];  			// Adimensional size of ROIs
		double[][][] Xnum    	= new double[(int)maxsize+1][q.length][3]; 	// Xnum(q)
		double deltavalue 	= 0.01;
		double [] delta 	= {0,-deltavalue,deltavalue};

		double P = 0;								// Acum probability
		

		GenericDialog gd = new GenericDialog("Select min and max size of box (in pixels)");	// Create dialog
		gd.addNumericField("minsize: ", 1, 0);							// Field1
		gd.addNumericField("maxsize: ", maxsize, 0);						// Field2
		gd.showDialog();									// Show
		if (gd.wasCanceled()) return;								// If no input, cancell program

		int nmin       = (int)gd.getNextNumber();						// Save value1
		int nmax       = (int)gd.getNextNumber();						// Save value2

		Xnum[nmax][1][1]=0;
		// Box Counting bucle
		// n=nmin; n<=nmax; n=n+A donde A puede ser, 1,2,4,8,16... en pedometrics cogi 16 
		for (int n = nmin; n<=nmax; n = n + 16){					// for each size
			IJ.log(""+(nmax-n));
			epsilon[n]   = (double)(n)/(double)maxsize;			// Adimensional size
			int numpixel = (n);						// Size in pixels of ROIs
			int box      = 0;
			double boxnumbers  = (woriginal-n+1)*(horiginal-n+1)*(doriginal-n+1);	// Number of box for gliding
			double[] Parchive  = new double[(int)boxnumbers];


			// Creating ROIs along image for that size
			for (int z=1; z<=stack.getSize()-numpixel+1; z++){
				IJ.log("sub"+(stack.getSize()-numpixel+1-z));
				for (int y=fullimage.y; y<fullimage.x+fullimage.width-numpixel+1; y++){
					for (int x=fullimage.x; x<fullimage.x+fullimage.width-numpixel+1; x++){
						P = 0;					// Inicialize P

						// Counting pixels in each slice of the box
						for (int zz=z; zz<=z+numpixel-1; zz++){
							salida_ip = stack.getProcessor(zz);	
							cuadradointeres = new Roi(x, y, numpixel, numpixel); 	// Size of ROI
							salida_ip.setRoi(cuadradointeres);			// Fix ROI in  image
							roi = salida_ip.getRoi();				// Processing ROI

							// Counting pixels along roi 
							for (int yy=roi.y; yy<(roi.y+roi.height); yy++){
								for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
									P += salida_ip.get(xx,yy)/sum;
								}
							}
						
						}
						// Processing dates of the box
						if (P>0) {
							Parchive[box] = P;
							box++;
							// Calculate Xnum
							for (int qposition = 0; qposition<q.length; qposition++){
								for (int deltaposition = 0; deltaposition<3; deltaposition++){
									Xnum[n][qposition][deltaposition] += pow(P,q[qposition]+delta[deltaposition]);
								}
							}	
						}						
					}
				}
			}
			for (int qposition = 0; qposition<q.length; qposition++){
				for (int deltaposition = 0; deltaposition<3; deltaposition++){
					Xnum[n][qposition][deltaposition] = Xnum[n][qposition][deltaposition]/box;
				}
			}
		}
		
		if (resolution == 16){
			Rectangle r = salida_ip.getRoi();
			for (int z=1; z<stack.getSize()+1; z++){
				salida_ip = stack.getProcessor(z);
      				for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){			// Normalize matriz
					for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
          					salida_ip.set(x, y, (int)abs(65535 - salida_ip.get(x,y)));
    					}	
  				}
			}
		imp.updateAndDraw();
		IJ.wait(500); 
		}
		else if (resolution != 16){
			IJ.run("Invert", "stack");				// Invert image. Count black mass
			IJ.wait(500);		
		}

		double minx0 = 0;
		double maxx0 = 0;
		double[] x = new double[nmax-nmin+1];				// x = log(epsilon)
		for(int n = 0; n <= nmax-nmin; n++){
			x[n] = log(epsilon[n+nmin]);
			//if (x[n]>maxx0){maxx0=x[n];}
			//if (x[n]<minx0){minx0=x[n];}
		}

		double miny0 = 0;
		double maxy0 = 0;	
		double[][][] y = new double[nmax-nmin+1][q.length][3]; 				// y = log(X(q))
		for (int n = 0; n<=nmax-nmin; n++){
			for (int qposition = 0; qposition<q.length; qposition++){
				for (int deltaposition = 0; deltaposition<3; deltaposition++){
					y[n][qposition][deltaposition] = log(Xnum[n+nmin][qposition][deltaposition]);
					//if (y[n][qposition][deltaposition]>maxy0){maxy0=y[n][qposition][deltaposition];}
					//else if (y[n][qposition][deltaposition]<miny0){miny0=y[n][qposition][deltaposition];}
				}
			}
		}

		// Intetanr limpiar valores de x e y
		int newSize = 0;
		for (int i = 0; i < x.length; i++)    {
		    if (!Double.isInfinite(x[i]) && !Double.isNaN(x[i])) {
			newSize++;
		    }
		}

		// Populate the new array.
		double[] x_new = new double[newSize];
		int newIndex = 0;
		for (int i = 0; i < x.length; i++) {
		    if (!Double.isInfinite(x[i])) {
			x_new[newIndex] = x[i];
			newIndex++;
		    }
		}

		// Limpiar y es mas complicado
		newSize = 0;
		for (int i = 0; i < x.length; i++)    {
		    if (!Double.isInfinite(y[i][0][0]) && !Double.isNaN(y[i][0][0])) {
			newSize++;
		    }
		}

		// Populate the new array.
		double[][][] y_new  = new double[newSize][q.length][3]; 	
		newIndex = 0;
		for (int i = 0; i < x.length; i++) {
			if (!Double.isInfinite(y[i][0][0]) && !Double.isNaN(y[i][0][0])) {
				for (int qposition = 0; qposition<q.length; qposition++){
					for (int deltaposition = 0; deltaposition<3; deltaposition++){
						y_new[newIndex][qposition][deltaposition] = y[i][qposition][deltaposition];
						if (y_new[newIndex][qposition][deltaposition]>maxy0){maxy0=y_new[newIndex][qposition][deltaposition];}
						if (y_new[newIndex][qposition][deltaposition]<miny0){miny0=y_new[newIndex][qposition][deltaposition];}
		    			}
				}
				newIndex++;
			}
		}


		for(int i = 0; i < x_new.length; i++){
			if (x_new[i]>maxx0){maxx0=x_new[i];}
			if (x_new[i]<minx0){minx0=x_new[i];}
		}


		// Return the new array.
		
		// CURVE FITTING (type of fitType: STRAIGHT_LINE=0, POLY2=1, POLY3=2, POLY4=3, EXPONENTIAL=4, POWER=5,...)	
		double[] TauQ     	= new double[q.length];					// Slopes
		double[][] TauQaux     	= new double[q.length][3];				// Slopes
		double[] desvTauQ 	= new double[q.length];					// Tau error
		double[][] desvTauQaux 	= new double[q.length][3];
		double[] relerror	= new double[q.length];
		double[] Dq       	= new double[q.length];					// Dq
		double[] desvDq   	= new double[q.length];					// Dq error
		double[] alphaq   	= new double[q.length];					// Dq error
		double[] desvalphaq   	= new double[q.length];
		double[] fq  		= new double[q.length];					// Dq error
		double[] desvfq  	= new double[q.length];
		double[] ybar 		= new double[y.length];
		
		// TauQ and eror
		for(int qposition = 0; qposition < q.length; qposition++){
			double[] yaux = new double[y_new.length];					// Adjusting aux array
			for (int deltaposition = 0; deltaposition<3; deltaposition++){
				for(int n = 0; n < y_new.length; n++){
					yaux[n] = y_new[n][qposition][deltaposition];
				}
				linefit fit1 = new linefit();
				fit1.x   = x_new;							// slope parameter x
				fit1.y   = yaux;						// slope parameter y
				TauQaux[qposition][deltaposition]     = fit1.getslope() - 3;	// slope
				desvTauQaux[qposition][deltaposition] = fit1.getdevslope();	// slope error
				if (deltaposition == 0){
			  		TauQ[qposition]     	= fit1.getslope() - 3;
					desvTauQ[qposition]	= fit1.getdevslope();		// slope error
				}
			}
			relerror[qposition] = desvTauQ[qposition]/abs(TauQ[qposition]);
			if (q[qposition] == 1){
				relerror[qposition] = (relerror[qposition-1]+relerror[qposition+1])/2;
			}	  	
		}

		// Dq and error
		for (int qposition = 0; qposition<q.length; qposition++){
			if (q[qposition] != 1){
				Dq[qposition]     = TauQ[qposition]/(q[qposition] - 1);
				desvDq[qposition] = desvTauQ[qposition]/abs(q[qposition] - 1);
			}	
		}
		for (int qposition = 0; qposition<q.length; qposition++){
			if (q[qposition] == 1){
				Dq[qposition]     = (Dq[qposition-1]+Dq[qposition+1])/2;
				desvDq[qposition] = (desvDq[qposition-1]+desvDq[qposition+1])/2;
			}
		}
		
		new WaitForUserDialog("CCCC", "Presionar para mostrar las figuras").show();
		IJ.run("Duplicate...","title=Los resultados corresponden a esta.tif duplicate");	// Duplicate image		
		
		// PLOT AND TABLES //
		String str = new String(" 3D Multifractal Gliding method");

		// PLot X(q,epsilon) vs epsilon
		Plot plot0 = new Plot("X(q,epsilon)"+str, "ln(epsilon)", "ln(X(q))");
		PlotWindow.noGridLines = false; 		// draw grid lines
		plot0.setLimits(minx0-0.5, maxx0 + 0.5,miny0-1,maxy0+1);
		for(int qposition = 0; qposition < q.length; qposition++){
			for(int n = 0; n < y_new.length; n++){
				ybar[n] = y_new[n][qposition][0];
			}
		plot0.addPoints(x_new, ybar, Plot.LINE);
		}
		plot0.show();
		ImagePlus imp_plot0 = plot0.getImagePlus();
		IJ.saveAs(imp_plot0, "Jpeg", path + "X(q,epsilon)"+ str);
		imp_plot0.close();
		IJ.selectWindow("X(q,epsilon)"+ str);	

		// Plot TauQ vs q
		PlotWindow.noGridLines = false; 						// draw grid lines
		Plot plot1 = new Plot("TauQ vs q"+str, "q", "TauQ", q, TauQ);
		plot1.setLimits((q[0]) - 0.5, q[q.length -1] + 0.5, TauQ[0]-desvTauQ[0] - 0.5, TauQ[q.length -1] + desvTauQ[0] + 0.5);
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
		plot2.setLimits((q[0]) - 0.5, q[q.length -1] + 0.5, Dq[q.length -1] - desvDq[q.length -1] - 0.5, Dq[0] + desvDq[0] + 0.5);
		plot2.setLineWidth(2);								// Line width
		plot2.addErrorBars(desvDq); 							// Errors bars
		plot2.setColor(Color.red);							// line color
		plot2.show();
		ImagePlus imp_plot2 = plot2.getImagePlus();
		IJ.saveAs(imp_plot2, "Jpeg", path + "Dq vs q"+ str);
		imp_plot2.close();

		// We estimate alpha(q) with central differences methods. 
		// alpha(q) = dTau(q)/dq;
		double maxalpha = 0;
		double minalpha = 10;
		double maxdesvalphaq 	= 0;

		for (int qposition = 0; qposition<q.length; qposition++){
			alphaq[qposition] = (TauQaux[qposition][2]-TauQaux[qposition][1])/((q[qposition]+deltavalue) - (q[qposition]-deltavalue));
			if (alphaq[qposition]<minalpha){minalpha=alphaq[qposition];}
			if (alphaq[qposition]>maxalpha){maxalpha=alphaq[qposition];}
			desvalphaq[qposition] = abs((alphaq[qposition])*(relerror[qposition]));
			if (desvalphaq[qposition]>maxdesvalphaq){maxdesvalphaq=desvalphaq[qposition];}
		}

		// Plot q vs alpha(q)
		PlotWindow.noGridLines = false; 						// draw grid lines
		Plot plot3 = new Plot("alpha(q) vs q"+str, "q", "alpha(q)", q, alphaq);
		plot3.setLimits(q[0] -0.5, q[q.length -1] + 0.5, alphaq[q.length -1] - 0.5 - maxdesvalphaq, alphaq[0] + 0.5 +  maxdesvalphaq );
		plot3.setLineWidth(2);								// Line width
		plot3.addErrorBars(desvalphaq); 							// Errors bars
		plot3.setColor(Color.red);							// line color
		//plot3.show();

		// Calculate f(alpha)
		// f(alpha) = alpha(q)*q - Tauq(q)
		double maxy = 0;
		double miny = 10;
		double maxdesvfq 	= 0;

		for (int qposition = 0; qposition<q.length; qposition++){
			fq[qposition] = alphaq[qposition]*q[qposition] - TauQ[qposition];
			if (fq[qposition] < miny){miny = fq[qposition];}
			if (fq[qposition] > maxy){maxy = fq[qposition];}
			desvfq[qposition] = (desvalphaq[qposition])*(abs(q[qposition])) + desvTauQ[qposition];
			if (desvfq[qposition]>maxdesvfq){maxdesvfq=desvfq[qposition];}
		}	  
		
		// Plot alpha vs f(alpha)
		PlotWindow.noGridLines = false; 						// draw grid lines
		Plot plot4 = new Plot("f(alpha) vs alpha"+str, "alpha", "f(alpha)", alphaq, fq);
		plot4.setLimits(minalpha - 0.5 - maxdesvalphaq,maxalpha + 0.5 + maxdesvalphaq,miny -0.5 - maxdesvfq,maxy + 0.5 + maxdesvfq);
		plot4.setLineWidth(2);								// Line width
		plot4.setColor(Color.red);							// line color
		for (int qposition = 0; qposition<q.length; qposition++){
			plot4.drawLine(alphaq[qposition]-desvalphaq[qposition],fq[qposition],alphaq[qposition]+desvalphaq[qposition],fq[qposition]);
		}
		plot4.show();
		ImagePlus imp_plot3 = plot4.getImagePlus();
		IJ.saveAs(imp_plot3, "Jpeg", path + "f(alpha) vs alpha"+ str);
		imp_plot3.close();

		// RESULTS TABLE
 		ResultsTable table = new ResultsTable();
 		for(int counter = 0; counter < q.length; counter++)    {
			table.incrementCounter();
			table.addValue("q", IJ.d2s(q[counter],2));
			table.addValue("TauQ",IJ.d2s(TauQ[counter], 8));
			table.addValue("error TauQ",IJ.d2s(desvTauQ[counter],8));
			table.addValue("Dq", IJ.d2s(Dq[counter],8));
			table.addValue("error Dq", IJ.d2s(desvDq[counter], 8));
 			table.addValue("alpha", IJ.d2s(alphaq[counter],8));
 			table.addValue("error alpha", IJ.d2s(desvalphaq[counter],8));
 			table.addValue("f(alpha)", IJ.d2s(fq[counter],8));
			table.addValue("error f", IJ.d2s(desvfq[counter],8));
		}
		table.updateResults();
		table.show("Results");
		IJ.selectWindow("Results");
		IJ.save(path+"Results.txt");

		time_end = System.currentTimeMillis();
		IJ.log("");
		IJ.log("The task has taken "+ ( time_end - time_start ) +" milliseconds");

	}

}
