// Ivan Gonzalez Torre
// This plugin realizes the monofractal analysis of a gray image
// It selects square region of interest from the original image
// First it converts the gray scale image to a Black and White one. 
// Here we implements the parameter which will discriminalize for doing the byte analysis
// Then calculates the monofractal Box Counting
// Finally it shows the fractal dimension and a plot with the dates and fitting curve.

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
import static java.lang.Math.*;					// Using static it is not necesary to call the class


public class _2D_Box_Counting_Monofractal implements PlugInFilter {

	public int setup(String arg, ImagePlus imp) {
		if (imp.getProcessor().isInvertedLut()){
			IJ.run("Invert LUT");
		}
		return DOES_8G+DOES_16+DOES_32;			// Suports gray scale image
	}

	public void run(ImageProcessor ip) {			// Process image and save in IP
		
		// CHECK SIZE OF IMAGE AND RESIZE IF NECESARY
		double woriginal = ip.getWidth();		// width from the original image
		double w = floor(log(woriginal)/log(2));	// log2 from the width
		double horiginal  = ip.getHeight();		// heigth from the original image
		double h = floor(log(horiginal)/log(2));	// log2 from the heigth
		w = h = min(pow(2,h),pow(2,w));			// new sizes of the image (it is square)
		int choose;
		ImagePlus salida = NewImage.createByteImage("Analyzed image", (int)w, (int)h, 1, NewImage.FILL_BLACK); 
		if (woriginal!=w || horiginal!=h){
			GenericDialog gd = new GenericDialog("Cut or resize image");				// Create dialog
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
				Roi subpart = new Roi(0, 0, w, h);
				ip.setRoi(subpart);
				IJ.run("Copy", "");
				salida = NewImage.createByteImage("Analyzed image", (int)w, (int)h, 1, NewImage.FILL_BLACK);
				IJ.run("Options...", "iterations=1 black count=1");
				salida.show();									// Show-
				IJ.run("Paste", "");								// Copy pixels
				salida.show();									// Show-
			}
			else if (choose == 1){	
				// Resize image
				salida = NewImage.createByteImage("Analyzed image", (int)woriginal, (int)horiginal, 1, NewImage.FILL_BLACK);
				ImageProcessor salida_aux = salida.getProcessor(); 				// Processing
				salida_aux.copyBits(ip,0,0,Blitter.COPY);					// Copy pixels	
				salida.show();									// Show
				salida.updateAndDraw();								// Update
				IJ.run("Size...", "width="+w+" height="+h+" average interpolation=Bilinear");	// Resize image
			}
		}
		else if (woriginal==w && horiginal==h) {
			salida = NewImage.createByteImage("Analyzed image", (int)w, (int)h, 1, NewImage.FILL_BLACK);
			ImageProcessor salida_aux = salida.getProcessor(); 					// Processing
			salida_aux.copyBits(ip,0,0,Blitter.COPY);						// Copy pixels	
			salida.show();										// Show
			salida.updateAndDraw();									// Update
		}

		// Threshold for converting the image to B&W
		IJ.run("Convert to Mask");
			
		// Dialog for continue
		GenericDialog gd1	= new GenericDialog("Process image");				// Create dialog
		String[] titleArray2 	= new String[2];
		gd1.addMessage("Occupied boxes wil be the black ones. Do you want to analyze this image or the inverted one?");
		titleArray2[0] = "This image";
		titleArray2[1] = "Inverted image";
		gd1.addChoice("  ", titleArray2, titleArray2[0]);
		gd1.showDialog();
		if (gd1.wasCanceled()) return;
		choose = gd1.getNextChoiceIndex();
		if (choose == 0){
		}
		else if (choose == 1){
			IJ.run("Invert");
			IJ.wait(500);
		}
		ImageProcessor salida_ip = salida.getProcessor(); 		// Processing the image after the changes
		
		long time_start, time_end;
		time_start = System.currentTimeMillis();

		// Begining of the program
		int epsilonmax      =  (int)(floor(log(w)/log(2)));		// Max times the image can be divided
		Rectangle fullimage = salida_ip.getRoi();			// Dates from the image
		Roi cuadradointeres;						// variable tipo Roi
		Rectangle roi;							// Region of interest multifractal
		
		GenericDialog gd = new GenericDialog("Select min and max size of box (in pixels)");	// Create dialog
		gd.addNumericField("minsize: ", 1, 0);							// Field1
		gd.addNumericField("maxsize: ", pow(2,epsilonmax), 0);					// Field2
		gd.showDialog();									// Show
		if (gd.wasCanceled()) return;		
		int pixelintromin       = (int)gd.getNextNumber();					// Save value1
		int pixelintromax       = (int)gd.getNextNumber();					// Save value2

		int pixelusedmin = (int)pow(2,(int)(log(pixelintromin)/log(2)));
		int pixelusedmax = (int)pow(2,(int)(log(pixelintromax)/log(2)));

		IJ.log("The minimum size of box that will be used is: "+ pixelusedmin);
		IJ.log(" ");
		IJ.log("The maximum size of box that will be used is: "+ pixelusedmax);

		epsilonmax = (int)((log(fullimage.width/pixelusedmin))/log(2));
		int epsilonmin = (int)((log(fullimage.width/pixelusedmax))/log(2));

		double[] epsilon    = new double[epsilonmax+1-epsilonmin];  			// Adimensional size of ROIs
		double[] N          = new double[epsilonmax+1-epsilonmin];			// Number of black ROI	

		// Box Counting bucle
		for (int n = 0; n<=epsilonmax-epsilonmin; n++){
			epsilon[n]   =  1/pow(2,n+epsilonmin);		
			int numpixel = fullimage.width/(int)pow(2,n+epsilonmin);		//Size in pixels of ROIs
			
			// Creating ROIs along image
			for (int y=fullimage.y; y<fullimage.y+fullimage.height; y=y+numpixel){
				for (int x=fullimage.x; x<fullimage.x+fullimage.width; x=x+numpixel){
					cuadradointeres = new Roi(x, y, numpixel, numpixel);	// Creating size of ROI
					salida_ip.setRoi(cuadradointeres);		    	// Fix ROI in the image
					roi = salida_ip.getRoi();				// Processing ROI
				
					// Check black pixels in ROI
					bucle1:
					for (int yy=roi.y; yy<(roi.y+roi.height); yy++){		
						for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
							if (salida_ip.get(xx,yy) == 0){		// If Black, count and exit
								N[n]++; 
								break bucle1;
							}
						}
					}
				}
			}
		}

		boolean isstack = false;
              	resultmonofractal result = new resultmonofractal();
		String str = new String(" 2D Monofractal Box Counting method");
		result.str = str;
		result.isstack = isstack;
		result.epsilon = epsilon;
		result.N = N;
		result.out();
		time_end = System.currentTimeMillis();
		IJ.log("");
		IJ.log("The task has taken "+ ( time_end - time_start ) +" milliseconds");
	}

}

