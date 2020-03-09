// This method calculates the configuration entropy function.
// Black and white image is required.

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

public class gray_lacunarity{
	ImageProcessor ip;
	String str;
	String path;

	public void out() {

		int woriginal  = ip.getWidth();				// width from the original image						
		int horiginal  = ip.getHeight(); 			// heigth from the original image
		int lmax = min(woriginal,horiginal);			// maximo cuadrado posible
		double[] lsize = new double[(int)lmax];
		double[] lacunarity = new double[(int)lmax];
		Roi cuadradointeres;
		Rectangle roi;				
		Rectangle fullimage = ip.getRoi();

		// Gliding Box bucle
		int lsize_iter = 0;
		for (int n = 1; n<=lmax; n++){ 					// Para cada tamano l de subcuadrado
			IJ.log(""+(lmax-n));					// Show bucles left 	
			int counter = 0;					// Contador de slices			
			int numpixel = (n);					// Size in pixels of ROIs
			int box      = 0;
			double boxnumbers  = (woriginal-n+1)*(horiginal-n+1);	// Number of posible boxes for gliding with box size n
			double[] measure_box  = new double[(int)boxnumbers];	// Vector donde se guarda el numero de boxes con value = 0 en cada slicing
			double[] measure_box_counting  = new double[(int)boxnumbers];	// Vector donde se guarda el numero de boxes con value = 0 en cada slicing															
			double numblackpixel;					// Contador de black boxes en cada slicing	
			double total_pixels=0;					

			// GLIDING FOR CALCULATING LACUNARITY
			for (int y=fullimage.y; y<fullimage.y+fullimage.height-numpixel+1; y++){
				for (int x=fullimage.x; x<fullimage.x+fullimage.width-numpixel+1; x++){	
					cuadradointeres = new Roi(x, y, numpixel, numpixel);	// Creating size of ROI
					ip.setRoi(cuadradointeres);		// Fix ROI in the image
					roi = ip.getRoi();			// Processing ROI
					int mass=0;					// Reinicia el contador de ese slice
					total_pixels = 0;				// Reinicia total pixels
					// Counting pixels along roi
					for (int yy=roi.y; yy<(roi.y+roi.height); yy++){		
						for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
							total_pixels += 255;
							mass += ip.get(xx,yy);
						}
					}
					
					measure_box[counter] = mass;	// Guarda el numero de pixels BLANCOS
					counter += 1;

				}
			}

			// Probability mass vector for each size l
			double[] freq_distribution = new double[(int)total_pixels + 1];	// vector p(black_boxes). Si el tamano de un cuadrado es 2x2, el posible numero de black boxes es 255*4+1 :(0,1,2,3,4)
			for (int pos = 0; pos<measure_box.length; pos++){		// Simplemente guarda el numero de boxes.
				freq_distribution[(int)measure_box[pos]] += 1;
			}


			for (int dum = 0; dum<freq_distribution.length; dum++){ 	// Convert in probabilities: Probability of box of size L with m black pixels
				freq_distribution[dum] = (double)freq_distribution[dum]/(double)boxnumbers;
			}
			
			// CALCULATE LACUNARITY AND CONFIGURATION ENTROPY
			double z1 = 0;
			double z2 = 0;		
			for (int m = 0; m<freq_distribution.length; m++){ 		// Iterate over the vector of probabilities
				if (freq_distribution[m] > 0){ 				// Para evitar log 0
					z1 += pow(m,1) * (double)freq_distribution[m];
					z2 += pow(m,2) * (double)freq_distribution[m];
				}
			}
			lacunarity[lsize_iter] = (double)z2/(pow(z1,2)); 		// Lacunarity vector	
			lsize[lsize_iter] = n;
			lsize_iter += 1;	

		}
		// RESULTS 
		String str = new String(" 2D Multifractal Gliding method");
		// LACUNARITY PLOT
		PlotWindow.noGridLines = false; 							// draw grid lines
		Plot plot_lacunarity = new Plot("lacunarity(l)", "l", "Lacunarity", lsize, lacunarity);
		plot_lacunarity.setLineWidth(2);							// Line width
		plot_lacunarity.setColor(Color.red);							// line color
		plot_lacunarity.show();

		ImagePlus imp_plot0 = plot_lacunarity.getImagePlus();
		IJ.saveAs(imp_plot0, "Jpeg", path + "Lacunarity"+ str);
		imp_plot0.close();


		// CURVE FITTING LACUNARITY (type of fitType: STRAIGHT_LINE=0)	
		double[] x = new double[90];					// x = log(1/epsilon)
		double[] y = new double[90];						// y = log(N)
		for(int counter = 0; counter < 90; counter++){
                	x[counter] = log(lsize[counter+10]);
			y[counter] = log(lacunarity[counter+10]);
                }

		linefit fit1 = new linefit();
		double slope, devslope;
		CurveFitter cf;
		fit1.x   = x;									// slope parameter x
		fit1.y   = y;									// slope parameter y
		slope    = fit1.getslope();							// slope
		devslope = fit1.getdevslope();							// slope error
		

		// RESULTS TABLE
 		ResultsTable table = new ResultsTable();
 		for(int dumm_c = 0; dumm_c < lsize.length; dumm_c++)    {
			table.incrementCounter();
			table.addValue("lsize", IJ.d2s(lsize[dumm_c], 2));
			table.addValue("lacunarity",IJ.d2s(lacunarity[dumm_c], 8));
		}

		table.setPrecision(8);
		table.updateResults();
		table.setPrecision(8);

		table.show("Results");
		IJ.selectWindow("Results");
		IJ.save(path+"gray_lacunarity_values.txt");
		IJ.run("Close");

		
 		ResultsTable table2 = new ResultsTable();
		table2.incrementCounter();
		table2.addValue("Lacunarity", IJ.d2s(slope, 5));
		table2.addValue("Lacunarity error", IJ.d2s(devslope, 5));

		table2.setPrecision(8);
		table2.updateResults();
		table2.setPrecision(8);

		table2.show("Results2");
		IJ.selectWindow("Results2");
		IJ.save(path+"Lacunarity.txt");
		IJ.run("Close");
		IJ.selectWindow("Results");
		IJ.run("Close");


	}
}
