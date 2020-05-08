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

import ij.plugin.filter.PlugInFilter;
import java.awt.*;
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.io.File;
import static java.lang.Math.*;

/**
 * Tool that computes multifractal analysis using Box Counting sliding method
 * on 2D Gray images
 * 
 * @author Ivan G Torre
 */

public class _2DMultifrBC implements PlugInFilter {
    public ImagePlus salida;
    public String imagename;    
    public String base_path;
    public int setup(String arg, ImagePlus imp) {
        if (imp.getProcessor().isInvertedLut()){
            IJ.run("Invert LUT");
        }
        base_path = imp.getOriginalFileInfo().directory;
        imagename = imp.getOriginalFileInfo().fileName;
        return DOES_8G+DOES_16+DOES_32;
        
    }


    public void run(ImageProcessor ip) {
		IJ.log(" 	MULTIFRAC  Copyright (C) <2020>  <ivangtorre>\n" +
				"	License GNU General Public License 3\n"+
				"	This program comes with ABSOLUTELY NO WARRANTY\n" + 
				"   This is free software, and you are welcome to redistribute it\n" + 
				"   under certain conditions. Please cite:\n" +
				"	I. G. Torre and A. M. Tarquis, MULTIFRAC: An ImageJ plugin for multifractal and \n" + 
				"	multiscale characterization of 2D and 3D stack images");
    	
    	/** ----------Dialogs and options------------- */         
    	Qdialog ask_q = new Qdialog();
        double[] q = ask_q.getq();

        GenericDialog gbase = new GenericDialog("Select the base of the squares");
        gbase.addNumericField("base: ", 2, 0); 
        gbase.showDialog(); 
        if (gbase.wasCanceled()) return;        
        int base;
        base = (int)gbase.getNextNumber();

        int pixelintro;
        int pixelused;
        GenericDialog gd2 = new GenericDialog("Select min size of box"); 
        gd2.addNumericField("Minimum size of box (in pixels: ", 1, 0);
        gd2.showDialog();
        if (gd2.wasCanceled()) return;      
        pixelintro       = (int)gd2.getNextNumber();  

    	/** ---------Check size and resize if necessary------------- */         
        double woriginal = ip.getWidth();                
        double wexact = log(woriginal)/log(base);  
        double w = round(log(woriginal)/log(base));
        if (abs(wexact - w) < 0.000001) { w = round(w); } else {w = floor(log(woriginal)/log(base));}
        double horiginal  = ip.getHeight();               
        double hexact = log(horiginal)/log(base);    
        double h = round(log(horiginal)/log(base));
        if (abs(hexact - h) < 0.000001) { h = round(h); } else {h = floor(log(horiginal)/log(base));}
        
        w = h = min(pow(base,h),pow(base,w));
        
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


    	/** ---------Choose which color is max value------------- */         
        GenericDialog gd1   = new GenericDialog("Process image");
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
                //Rectangle r = salida_ip.getRoi();
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
            //Rectangle r = salida_ip.getRoi();
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
        
    	/** ---------Handle output path---------------------------------- */         
        String save_path_pred;
        String[] parts = imagename.split("\\.");
        String part1 = parts[0]; // 004
        base_path = base_path + part1 + "/";
        if (choose == 0){save_path_pred = base_path + "box_counting/white_max_value/";} 
        else {save_path_pred = base_path + "box_counting/black_max_value/";}
        String path;
        path = save_path_pred;

    	/** ---------Check if path exists and if not, create it----------------------------- */         
        File wdir = new File(path);
        if (!wdir.isDirectory() || !wdir.exists() || !wdir.canRead()) {
            wdir.mkdirs();
            IJ.log("Path created.");}


        long time_start, time_end;
        time_start = System.currentTimeMillis();

        
    	/** ---------Start of computational analysis---------------------------------- */         
        double sum = salida_ip.getStatistics().area*salida_ip.getStatistics().mean;     // Sum of all pixels value

        double[][] matriz = new double[fullimage.height][fullimage.width];      // Matriz with values
        for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){         // Normalize matriz
            for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
                matriz[x][y] = salida_ip.get(x,y)/sum;
            }
        }
        

        pixelused = (int)pow(base,(int)(log(pixelintro)/log(base)));
        IJ.log("The minimum size of box that will be used is: "+ pixelused);

        int epsilonmax = (int)((log(fullimage.width/pixelused))/log(base));
        double epsilonexact = (log(fullimage.width/pixelused))/log(base);
        if (abs(epsilonexact - round(epsilonexact)) < 0.000001) {epsilonmax = (int)(round(epsilonexact)); }

        IJ.log("epsilonmax: "+ epsilonmax);
        double[] epsilon     = new double[epsilonmax+1];               
        double[][] Xnum      = new double[epsilonmax+1][q.length];         
        double[][] denomMu       = new double[epsilonmax+1][q.length];
        double[][] numalpha      = new double[epsilonmax+1][q.length];
        double[][] numf          = new double[epsilonmax+1][q.length];

        double P = 0;                                  
        Roi cuadradointeres;                               
        Rectangle roi;                            
    
    	/** ---------Box Counting Loop---------------------------------- */         
        for (int n = 0; n<=epsilonmax; n++){                   
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
                    
                	/** ---------Xnum---------------------------------- */         
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
        	
            /** ---------Mu---------------------------------- */         
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
        //Rectangle r = salida_ip.getRoi();
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

    	/** ---------Handle Results---------------------------------- */         
		ToolMultifr result = new ToolMultifr();
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
		IJ.log("The task has taken "+ ( time_end - time_start ) +" milliseconds");

	}

}

