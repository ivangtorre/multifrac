// Ivan Gonzalez Torre

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
import static java.lang.Math.*;					// Using static it is not necesary to call the class


public class _3D_Box_Counting_Monofractal implements PlugInFilter {

protected ImageStack stackaux;

	public int setup(String arg, ImagePlus impaux) {
		stackaux = impaux.getStack();
		if (impaux.getProcessor().isInvertedLut()){
			IJ.run("Invert LUT");
		}
		
		return DOES_8G+DOES_16+DOES_32+STACK_REQUIRED;	// Suports gray scale image
	}


	public void run(ImageProcessor ip1) {			// Process image and save in IP
		
		//Resize image
		double woriginal = stackaux.getWidth();		// width from the original image
		double w = floor(log(woriginal)/log(2));	// log2 from the width
		double horiginal  = stackaux.getHeight();	// heigth from the original image
		double h = floor(log(horiginal)/log(2));	// log2 from the heigth
		double doriginal = stackaux.getSize();		// depth from the original image
		double d = floor(log(doriginal)/log(2));		
		w = h = min(h,w);				// new sizes of the image (it is square)
		w = h = d = min(pow(2,h),pow(2,d));

		// Creating the new image
		IJ.run("Options...", "iterations=1 black count=1"); 
		IJ.run("Duplicate...","title=Analyzed_image.tif duplicate");					// Duplicate image
		IJ.run("Size...", "width="+w+" height="+h+" depth="+d+" average interpolation=Bilinear");	// Resize image
		IJ.run("Options...", "iterations=1 black count=1");			// Check not invertign LUT	
		
		IJ.run("Convert to Mask","method=Default background=Dark black");	// Threshold for converting the image to B&W
		
		// Dialog for continue
		GenericDialog gd1	= new GenericDialog("Process image");				// Create dialog
		String[] titleArray2 	= new String[2];
		gd1.addMessage("Occupied boxes wil be the black ones. Do you want to analyze this image or the inverted one?");
		titleArray2[0] = "This image";
		titleArray2[1] = "Inverted image";
		gd1.addChoice("  ", titleArray2, titleArray2[0]);
		gd1.showDialog();
      		if (gd1.wasCanceled()) return;
		int choose;
		choose = gd1.getNextChoiceIndex();
		if (choose == 0){
		}
		else if (choose == 1){
			IJ.run("Invert", "stack");
			IJ.wait(500);
		}

		ImagePlus imp 		= IJ.getImage(); 
		ImageStack stack 	= imp.getStack(); 

		Roi cuadradointeres;							// variable tipo Roi
		Rectangle roi;								// Region of interest multifractal
		ImageProcessor salida_fullimage = stack.getProcessor(1);		// Processor first slice
		Rectangle fullimage = salida_fullimage.getRoi(); 			// characters of that slice
		ImageProcessor salida_ip;

		long time_start, time_end;
		time_start = System.currentTimeMillis();

		// Begining of the program
		int epsilonmax      = (int)(floor(log(w)/log(2)));			// Max times the image can be divided

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
		IJ.log(" ");

		epsilonmax = (int)((log(fullimage.width/pixelusedmin))/log(2));
		int epsilonmin = (int)((log(fullimage.width/pixelusedmax))/log(2));

		double[] epsilon    = new double[epsilonmax+1-epsilonmin];  		// Adimensional size of ROIs
		double[] N          = new double[epsilonmax+1-epsilonmin];		// Number of black ROI	

		// Box Counting bucle
		for (int n = 0; n<=epsilonmax-epsilonmin; n++){
			IJ.log(""+(epsilonmax-epsilonmin-n));
			epsilon[n]   =  1/pow(2,n+epsilonmin);		
			int numpixel = fullimage.width/(int)pow(2,n+epsilonmin);	//Size in pixels of ROIs
			
			// Creating 3D ROIs along image
			for (int z=1; z<stack.getSize()+1; z=z+numpixel){ 
				for (int y=fullimage.y; y<fullimage.y+fullimage.height; y=y+numpixel){
					for (int x=fullimage.x; x<fullimage.x+fullimage.width; x=x+numpixel){

						// Check black pixels in ROI
						bucle1:
						for (int zz=z; zz<z+numpixel; zz++){
							salida_ip = stack.getProcessor(zz);
							cuadradointeres = new Roi(x, y, numpixel, numpixel);	// Creating size of ROI
							salida_ip.setRoi(cuadradointeres);			// Fix ROI in the image
							roi = salida_ip.getRoi();				// Processing ROI
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
			}
		}
		IJ.log(""+epsilonmax);
		boolean isstack = true;
		resultmonofractal result = new resultmonofractal();
		String str = new String(" 3D Monofractal Box Counting method");
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

