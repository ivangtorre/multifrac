// Ivan G Torre
// ivan.gonzalez.torre@upm.es
// https://github.com/ivangtorre/multifractal-analysis-images
// Cite: Torre, I. G., Losada, J. C., Heck, R. J., & Tarquis, A. M. (2018). Multifractal analysis of 3D images of tillage soil. Geoderma, 311, 167-174
// *******************************
// *******************************

// This plugin realizes the  analysis of a 3D gray image
// Calculate multifractal analysis of 3D images by box counting including Dq, alpha, f(alpha)
// *******************************

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
import static java.lang.Math.*;
import java.io.File;
import java.util.Arrays;


public class _3D_Box_Counting_Multifractal implements PlugInFilter {

	protected ImageStack stackaux;
	public String imagename;	
	public String base_path;
	public int setup(String arg, ImagePlus imp) {
		stackaux = imp.getStack();
		if (imp.getProcessor().isInvertedLut()){
			IJ.run("Invert LUT");
		}
		base_path = imp.getOriginalFileInfo().directory;
		imagename = imp.getOriginalFileInfo().fileName;
		return DOES_8G+DOES_16+DOES_32+STACK_REQUIRED;
	}


	public void run(ImageProcessor ip) {			
		//INPUT DIALOG Q VALUES ********************************************************************************
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
			
		// CHECK SIZE OF IMAGE ***************************************************************************************
		double woriginal = ip.getWidth();				// Width from the original image
		double w = floor(log(woriginal)/log(2));			// log2 from the width
		double horiginal  = ip.getHeight();				// heigth from the original image
		double h = floor(log(horiginal)/log(2));			// log2 from the heigth
		double doriginal = stackaux.getSize();				// depth from the original image
		double d = floor(log(doriginal)/log(2));		
		w = h = min(h,w);						// new sizes of the image (it is square)
		w = h = d = min(pow(2,h),pow(2,d));


		// MAX AND MIN VALUES OF EPSILON ***************************************************************************
		// CHECK SIZE OF IMAGE AND RESIZE IF NECESARY ************************************************************
		int choose;
		int resolution = ip.getBitDepth();
		if (woriginal!=w || horiginal!=h || doriginal!=d){
			GenericDialog gd = new GenericDialog("Cut or resize image");	// Create dialog
			String[] titleArray2 = new String[2];
			gd.addMessage("The image has not the properly size, choose what to do");
			titleArray2[0] = "Cut image";
			titleArray2[1] = "Resize image";
			gd.addChoice("  ", titleArray2, titleArray2[0]);
			gd.showDialog();
      			if (gd.wasCanceled()) return;
			choose = gd.getNextChoiceIndex();
			if (choose == 0){
				// Cut image
				IJ.makeRectangle(0, 0, (int)w, (int)h);
				//Roi subpart = new Roi(0, 0, w, h);
				//ip.setRoi(subpart);
				IJ.run("Make Substack...", "slices=1-"+(int)d);
				IJ.run("Rename...", "title=Analyzed_image.tif");
			}
			else if (choose == 1){	
				// Resize image
				IJ.run("Options...", "iterations=1 black count=1"); 
				IJ.run("Duplicate...","title=Analyzed_image.tif duplicate");			// Duplicate image
				IJ.run("Size...", "width="+w+" height="+h+" depth="+d+" average interpolation=Bilinear");
				IJ.run("Options...", "iterations=1 black count=1");		// Check not invertign LUT
			}
		}
		else if (woriginal==w && horiginal==h) {
			// Copy image
			// Creating the new image
			IJ.run("Options...", "iterations=1 black count=1"); 
			IJ.run("Duplicate...","title=Analyzed_image.tif duplicate");			// Duplicate image
			IJ.run("Size...", "width="+w+" height="+h+" depth="+d+" average interpolation=Bilinear");// Resize image
			IJ.run("Options...", "iterations=1 black count=1");			// Check not invertign LUT
		}



		// OUTPUT PATH ****************************************************************************************************
		String save_path_pred;
		String[] parts = imagename.split("\\.");
		String part1 = parts[0]; // 004
		base_path = base_path + part1 + "/";
		//if (choose == 0){save_path_pred = base_path + "box_counting/white_max_value/";} 
		//else {save_path_pred = base_path + "box_counting/black_max_value/";}
		save_path_pred = base_path + "box_counting/";
		String path;
		path = save_path_pred;

		// Check if path exists and if not, create it
		File wdir = new File(path);
		if (!wdir.isDirectory() || !wdir.exists() || !wdir.canRead()) {
			wdir.mkdirs();
			IJ.log("Path created.");}



		// CONTINUE PROGRAM
		GenericDialog gd1	= new GenericDialog("Process image");			// Create dialog
		String[] titleArray2 	= new String[2];
		gd1.addMessage("Maximum measure will be for the black. Do you want to analyze this image or the inverted one?");
		titleArray2[0] = "This image";
		titleArray2[1] = "Inverted image";
		gd1.addChoice("  ", titleArray2, titleArray2[0]);
		gd1.showDialog();
      		if (gd1.wasCanceled()) return;
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
		int epsilonmax  =  (int)(floor(log(w)/log(2)));				// Max times the image can be divided

		// Choose the maximum and minimum values of epsilon
		int pixelintro;
		int pixelused;
		GenericDialog gd2 = new GenericDialog("Select min size of box");	// Create dialog
		gd2.addNumericField("Minimum size of box (in pixels: ", 1, 0);		// Field1
		gd2.showDialog();								// Show
      		if (gd2.wasCanceled()) return;
		pixelintro       = (int)gd2.getNextNumber();					// Save value1
		
		pixelused = (int)pow(2,(int)(log(pixelintro)/log(2)));
		IJ.log("The minimum size of box that will be used is: "+ pixelused);
		
		double sum = 0;								// Sum of all pixels value 
		for (int z=1; z<stack.getSize()+1; z++){ 
			ImageProcessor salidaux_ip = stack.getProcessor(z);
			sum += salidaux_ip.getStatistics().area*salidaux_ip.getStatistics().mean;
		}

		double P = 0;								// Acum probability

		epsilonmax = (int)((log(fullimage.width/pixelused))/log(2));
		double[] epsilon 	 = new double[epsilonmax+1];  			// Adimensional size of ROIs
		double[][] Xnum    	 = new double[epsilonmax+1][q.length]; 		// Xnum(q)
		double[][] denomMu       = new double[epsilonmax+1][q.length];
		double[][] numalpha    	 = new double[epsilonmax+1][q.length];
		double[][] numf     	 = new double[epsilonmax+1][q.length];

		// Box Counting bucle
		for (int n = 0; n<epsilon.length; n++){					// for each size
			IJ.log(""+(epsilon.length-n));
			epsilon[n]         = 1/pow(2,n);
			int numpixel 	   = fullimage.width/(int)pow(2,n);
			int box            = 0;
			double boxnumbers  = pow(2,n*2)*pow(2,n);
			double[] Parchive  = new double[(int)boxnumbers];


			// Creating ROIs along image for that size
			for (int z=1; z<stack.getSize()+1; z=z+numpixel){
				for (int y=fullimage.y; y<fullimage.y+fullimage.height; y=y+numpixel){
					for (int x=fullimage.x; x<fullimage.x+fullimage.width; x=x+numpixel){
						P = 0;					// Inicialize P	
						for (int zz=z; zz<z+numpixel; zz++){
							salida_ip = stack.getProcessor(zz);	
							cuadradointeres = new Roi(x, y, numpixel, numpixel); 	// Creating size of ROI
							salida_ip.setRoi(cuadradointeres);			// Fix ROI in  image
							roi = salida_ip.getRoi();				// Processing ROI

							// Counting pixels along roi 
							for (int yy=roi.y; yy<(roi.y+roi.height); yy++){
								for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
									P += salida_ip.get(xx,yy)/sum;
								}
							}
						
						}
						Parchive[box] = P;
						box++;

						// Calculate of Xnum
						for (int qposition = 0; qposition<q.length; qposition++){
							if (q[qposition] == 1){
								if (P>0){
									Xnum[n][qposition] += P*log(P);
								}
							}
							else{
								if (P>0){	
								Xnum[n][qposition] += pow(P,q[qposition]);
								}

							}
						}
						
										
					}
				}
			}
			// Calculate Mu
			for (int count = 0; count<box; count++){
				for (int qposition = 0; qposition<q.length; qposition++){
					if (Parchive[count]>0){denomMu[n][qposition] += pow(Parchive[count],q[qposition]);}
				}
			}
			for (int count = 0; count<box; count++){
				for (int qposition = 0; qposition<q.length; qposition++){
				if (Parchive[count]>0){numalpha[n][qposition] += (pow(Parchive[count],q[qposition])/denomMu[n][qposition])*log(Parchive[count]);
					
				numf[n][qposition] += (pow(Parchive[count],q[qposition])/denomMu[n][qposition])*log(pow(Parchive[count],q[qposition])/denomMu[n][qposition]);}
				}
			}
		}
		IJ.log("0");

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

		// Plot tables and result
              	resultmultifractal result = new resultmultifractal();
		String str = new String(" 3D Multifractal Box Counting method");
		result.str = str;
		boolean isstack = true;
		result.isstack = isstack;
		result.epsilon = epsilon;
		result.Xnum = Xnum;
		result.denomMu = denomMu;
		result.numalpha = numalpha;
		result.numf = numf;
		result.q = q;
		result.path = path;
		result.out();

		time_end = System.currentTimeMillis();
		IJ.log("");
		IJ.log("The task has taken "+ ( time_end - time_start ) +" milliseconds");
	}

}
