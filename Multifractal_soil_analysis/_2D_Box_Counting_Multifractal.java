// Ivan G Torre
// ivan.gonzalez.torre@upm.es
// https://github.com/ivangtorre/multifractal-analysis-images
// Cite: Torre, I. G., Martín-Sotoca, J. J., Losada, J. C., López, P., & Tarquis, A. M. (2020). Scaling properties of binary and greyscale images in the context of X-ray soil tomography. Geoderma, 365, 114205.
// *******************************
// *******************************

// This plugin realizes the  analysis of a 2D gray image
// Calculate multifractal analysis of 2D images by box counting method including Dq, alpha, f(alpha)
// *******************************

// Import libraries
import ij.plugin.filter.PlugInFilter;
import java.awt.*;
import ij.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.process.*;
import ij.text.*;
import ij.gui.*;
import ij.util.*;
import ij.io.*;
import ij.process.*;
import ij.measure.*;
import java.io.File;
import static java.lang.Math.*; 

public class _2D_Box_Counting_Multifractal implements PlugInFilter {
    public ImagePlus salida;
    public String imagename;    
    public String base_path;
    public int setup(String arg, ImagePlus imp) {
        if (imp.getProcessor().isInvertedLut()){
            IJ.run("Invert LUT");
        }
        base_path = imp.getOriginalFileInfo().directory;
        imagename = imp.getOriginalFileInfo().fileName;
        return DOES_8G+DOES_16+DOES_32;                     // Suports gray scale image
        
    }


    public void run(ImageProcessor ip) {                        // Process image and save in IP
        //INPUT DIALOG Q VALUES ********************************************************************************
         qdialog ask_q = new qdialog();
         double[] q = ask_q.getq();

        // INPUT DIALOG OF THE BASE *******************************************************************************
        GenericDialog gbase = new GenericDialog("Select the base of the squares");  // Create dialog
        gbase.addNumericField("base: ", 2, 0);                      // Field1
        gbase.showDialog();                             // Show
        if (gbase.wasCanceled()) return;        
        int base;
        base = (int)gbase.getNextNumber();                      // Save value1

        // MAX AND MIN VALUES OF EPSILON ***************************************************************************
        int pixelintro;
        int pixelused;
        GenericDialog gd2 = new GenericDialog("Select min size of box");        // Create dialog
        gd2.addNumericField("Minimum size of box (in pixels: ", 1, 0);          // Field1
        gd2.showDialog();                               // Show
        if (gd2.wasCanceled()) return;      
        pixelintro       = (int)gd2.getNextNumber();                    // Save value1

    
        // CHECK SIZE OF IMAGE AND RESIZE IF NECESARY ************************************************************
        double woriginal = ip.getWidth();                   // width from the original image
        double wexact = log(woriginal)/log(base);               // log2 from the width
        double w = round(log(woriginal)/log(base));
        if (abs(wexact - w) < 0.000001) { w = round(w); } else {w = floor(log(woriginal)/log(base));}
        double horiginal  = ip.getHeight();                 // heigth from the original image
        double hexact = log(horiginal)/log(base);               // log2 from the heigth
        double h = round(log(horiginal)/log(base));
        if (abs(hexact - h) < 0.000001) { h = round(h); } else {h = floor(log(horiginal)/log(base));}
        
        w = h = min(pow(base,h),pow(base,w));                   // new sizes of the image (it is square)
        
        int choose;
        int resolution = ip.getBitDepth();
        ImagePlus salida = NewImage.createImage("Analyzed Image", (int)w, (int)h, 1, resolution, NewImage.FILL_BLACK);
        if (woriginal!=w || horiginal!=h){
            GenericDialog gd = new GenericDialog("Cut or resize image");    // Create dialog
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
                salida = NewImage.createImage("Analyzed Image", (int)w, (int)h, 1, resolution, NewImage.FILL_BLACK);
                IJ.run("Options...", "iterations=1 black count=1");
                salida.show();                          // Show-
                IJ.run("Paste", "");                        // Copy pixels
                salida.show();                          // Show-
            }
            else if (choose == 1){  
                // Resize image
                salida = NewImage.createImage("Analyzed Image", (int)woriginal, (int)horiginal, 1, resolution, NewImage.FILL_BLACK);
                ImageProcessor salida_aux = salida.getProcessor();              // Processing
                salida_aux.copyBits(ip,0,0,Blitter.COPY);                   // Copy pixels  
                salida.show();                                  // Show
                salida.updateAndDraw();                             // Update
                IJ.run("Size...", "width="+w+" height="+h+" average interpolation=Bilinear");   // Resize image
            }
        }
        else if (woriginal==w && horiginal==h) {
            salida = NewImage.createImage("Analyzed Image", (int)woriginal, (int)horiginal, 1, resolution, NewImage.FILL_BLACK);
            ImageProcessor salida_aux = salida.getProcessor();                  // Processing
            salida_aux.copyBits(ip,0,0,Blitter.COPY);                       // Copy pixels  
            salida.show();                                      // Show
            salida.updateAndDraw();                         // Update
        }


        // CHOOSE THE WHICH PIXEL IS MAX VALUE ***************************************************************************
        GenericDialog gd1   = new GenericDialog("Process image");           // Create dialog
        String[] titleArray2    = new String[2];
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
            ImageProcessor salida_ip = salida.getProcessor();       // Processing the image after the changes
            Rectangle fullimage = salida_ip.getRoi();               // ROI full image
            if (resolution == 16){
                Rectangle r = salida_ip.getRoi();
                    for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){         // Normalize matriz
                    for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
                            salida_ip.set(x, y, (int)abs(65535 - salida_ip.get(x,y)));
                        }
                }
            salida.updateAndDraw();
            IJ.wait(500);
            }

            else if (resolution != 16){
                IJ.run(salida, "Invert", "");               // Invert image. Count black mass
                IJ.wait(500);
            }
        }

        ImageProcessor salida_ip = salida.getProcessor();           // Processing the image after the changes
        Rectangle fullimage = salida_ip.getRoi();               // ROI full image
        if (resolution == 16){
            Rectangle r = salida_ip.getRoi();
                for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){         // Normalize matriz
                for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
                    salida_ip.set(x, y, (int)abs(65535 - salida_ip.get(x,y)));
                    }
            }
        salida.updateAndDraw();
        IJ.wait(500);
        } 


        else if (resolution != 16){
            IJ.run(salida, "Invert", "");                       // Invert image. Count black mass
            IJ.wait(500);
        }
        
        // OUTPUT PATH ****************************************************************************************************
        String save_path_pred;
        String[] parts = imagename.split("\\.");
        String part1 = parts[0]; // 004
        base_path = base_path + part1 + "/";
        if (choose == 0){save_path_pred = base_path + "box_counting/white_max_value/";} 
        else {save_path_pred = base_path + "box_counting/black_max_value/";}

        String path;
        path = save_path_pred;

        // Check if path exists and if not, create it
        File wdir = new File(path);
        if (!wdir.isDirectory() || !wdir.exists() || !wdir.canRead()) {
            wdir.mkdirs();
            IJ.log("Path created.");}


        long time_start, time_end;
        time_start = System.currentTimeMillis();

        // Begining of the program
        double sum = salida_ip.getStatistics().area*salida_ip.getStatistics().mean;     // Sum of all pixels value

        double[][] matriz = new double[fullimage.height][fullimage.width];      // Matriz with values
        for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){         // Normalize matriz
            for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
                matriz[x][y] = salida_ip.get(x,y)/sum;
            }
        }
        

        pixelused = (int)pow(base,(int)(log(pixelintro)/log(base)));
        IJ.log("The minimum size of box that will be used is: "+ pixelused);

        //double epsilonexact = log(w)/log(base);               // log2 from the heigth
        //IJ.log("epsilonexact: "+ epsilonexact);
        //int epsilonmax = (int)(floor(log(w)/log(base)));
        //if (abs(epsilonexact - round(epsilonexact)) < 0.000001) {epsilonmax = (int)(round(epsilonexact)); } else {epsilonmax = (int)(floor(log(w)/log(base)));}

        int epsilonmax = (int)((log(fullimage.width/pixelused))/log(base));
        double epsilonexact = (log(fullimage.width/pixelused))/log(base);
        if (abs(epsilonexact - round(epsilonexact)) < 0.000001) {epsilonmax = (int)(round(epsilonexact)); } else {epsilonmax = epsilonmax;}

        IJ.log("epsilonmax: "+ epsilonmax);
        double[] epsilon     = new double[epsilonmax+1];                // Adimensional size of ROIs
        double[] N           = new double[epsilonmax+1];                // Number of black ROI
        double[][] Xnum      = new double[epsilonmax+1][q.length];          // Xnum(q)
        double[][] denomMu       = new double[epsilonmax+1][q.length];
        double[][] numalpha      = new double[epsilonmax+1][q.length];
        double[][] numf          = new double[epsilonmax+1][q.length];

        double P = 0;                                   // Acum probability
        Roi cuadradointeres;                                // variable tipo Roi
        Rectangle roi;                                  // Region of interest multifractal
    
        // Box Counting bucle
        for (int n = 0; n<=epsilonmax; n++){                    // for each size
            IJ.log("n: "+ n);           
            epsilon[n]         =  1/pow(base,n);
            int numpixel       = fullimage.width/(int)pow(base,n);
            int box            = 0;
            double boxnumbers  = pow(base,n*2)*pow(base,n);
            double[] Parchive  = new double[(int)boxnumbers];
            // Creating ROIs along image for that size
            for (int y=fullimage.y; y<fullimage.y+fullimage.height; y=y+numpixel){
                for (int x=fullimage.x; x<fullimage.x+fullimage.width; x=x+numpixel){
                    cuadradointeres = new Roi(x, y, numpixel, numpixel);    // Creating size of ROI
                    salida_ip.setRoi(cuadradointeres);              // Fix ROI in the image
                    roi = salida_ip.getRoi();               // Processing ROI
                    P = 0;                          // Inicialize P 
                    // Counting pixels along roi 
                    for (int yy=roi.y; yy<(roi.y+roi.height); yy++){        
                        for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
                            P += matriz[xx][yy];
                        }
                    }
                    Parchive[box] = P;
                    box++;
                    // Calculate of Xnum
                    for (int qposition = 0; qposition<q.length; qposition++){
                        if (q[qposition] == 1){
                            if(P>0){Xnum[n][qposition] += P*log(P);}
                        }
                        else{   
                            if(P>0){Xnum[n][qposition] += pow(P,q[qposition]);}
                        }
                    }               
                }
            }
            // Calculate Mu
            for (int count = 0; count<box; count++){
                for (int qposition = 0; qposition<q.length; qposition++){
                    if(Parchive[count]>0){denomMu[n][qposition] += pow(Parchive[count],q[qposition]);}
                }
            }
            for (int count = 0; count<box; count++){
                for (int qposition = 0; qposition<q.length; qposition++){
                if(Parchive[count]>0){numalpha[n][qposition] += (pow(Parchive[count],q[qposition])/denomMu[n][qposition])*log(Parchive[count]);
                    
                numf[n][qposition] += (pow(Parchive[count],q[qposition])/denomMu[n][qposition])*log(pow(Parchive[count],q[qposition])/denomMu[n][qposition]);}
                }
            }
             
        }
        if (resolution == 16){
        Rectangle r = salida_ip.getRoi();
                for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){         // Normalize matriz
                for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
                        salida_ip.set(x, y, (int)abs(65535 - salida_ip.get(x,y)));
                    }
            }
        salida.updateAndDraw();
        IJ.wait(500);       
        } 
        

		else if (resolution != 16){
			IJ.run(salida, "Invert", "");						// Invert image. Count black mass
			IJ.wait(500);		
		}

		// Plot tables and result
		resultmultifractal result = new resultmultifractal();
		String str = new String(" 2D Multifractal Box Counting method");
		result.str = str;
		boolean isstack = false;
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

