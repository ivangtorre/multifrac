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
import ij.gui.*;
import ij.process.*;
import static java.lang.Math.*;
import java.io.File;

/**
 * Tool that computes multifractal analysis using Box Counting sliding method
 * on 3D stack Gray images
 * 
 * @author Ivan G Torre
 */

public class _3DMultifrBC implements PlugInFilter {
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
		IJ.log(" 	MULTIFRAC  Copyright (C) <2020>  <ivangtorre>\n" +
				"	License GNU General Public License 3\n"+
				"	This program comes with ABSOLUTELY NO WARRANTY\n" + 
				"   This is free software, and you are welcome to redistribute it\n" + 
				"   under certain conditions. Please cite:\n" +
				"	Torre, I.G., Heck R.J. & Tarquis, A.M. (2020).\n" +
				"   MULTIFRAC: An ImageJ plugin for multiscale characterization\n" + 
				"   of 2D and 3D stack images.");
    	
    	/** ----------Dialogs and options------------- */   
        Qdialog ask_q = new Qdialog();
        double[] q = ask_q.getq();
            
    	/** ---------Check size and resize if necessary------------- */         
        double woriginal = ip.getWidth(); 
        double w = floor(log(woriginal)/log(2));
        double horiginal  = ip.getHeight();
        double h = floor(log(horiginal)/log(2));
        double doriginal = stackaux.getSize();
        double d = floor(log(doriginal)/log(2));        
        w = h = min(h,w);
        w = h = d = min(pow(2,h),pow(2,d));

        int choose;
        int resolution = ip.getBitDepth();
        if (woriginal!=w || horiginal!=h || doriginal!=d){
            GenericDialog gd = new GenericDialog("Cut or resize image");
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
                IJ.run("Make Substack...", "slices=1-"+(int)d);
                IJ.run("Rename...", "title=Analyzed_image.tif");
            }
            else if (choose == 1){  
                // Resize image
                IJ.run("Options...", "iterations=1 black count=1"); 
                IJ.run("Duplicate...","title=Analyzed_image.tif duplicate");
                IJ.run("Size...", "width="+w+" height="+h+" depth="+d+" average interpolation=Bilinear");
                IJ.run("Options...", "iterations=1 black count=1");
            }
        }
        else if (woriginal==w && horiginal==h) {
            IJ.run("Options...", "iterations=1 black count=1"); 
            IJ.run("Duplicate...","title=Analyzed_image.tif duplicate");
            IJ.run("Size...", "width="+w+" height="+h+" depth="+d+" average interpolation=Bilinear");
            IJ.run("Options...", "iterations=1 black count=1");
        }


    	/** ---------Handle output path----------------------------------------------------------- */         
        String save_path_pred;
        String[] parts = imagename.split("\\.");
        String part1 = parts[0];
        base_path = base_path + part1 + "/";
        save_path_pred = base_path + "box_counting/";
        String path;
        path = save_path_pred;

    	/** ---------Check if path exists and if not, create it----------------------------- */         
        File wdir = new File(path);
        if (!wdir.isDirectory() || !wdir.exists() || !wdir.canRead()) {
            wdir.mkdirs();
            IJ.log("Path created.");}


    	/** ---------Choose which color is max value------------- */         
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
            ImagePlus imp = IJ.getImage(); 
            ImageStack stack = imp.getStack(); 
            //Roi cuadradointeres;                            // variable tipo Roi
            //Rectangle roi;                              // Region of interest multifractal
            ImageProcessor salida_ip = stack.getProcessor(1);
            Rectangle fullimage = salida_ip.getRoi();

            if (resolution == 16){
                //Rectangle r = salida_ip.getRoi();
                for (int z=1; z<stack.getSize()+1; z++){
                    salida_ip = stack.getProcessor(z);
                        for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){
                        for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
                            salida_ip.set(x, y, (int)abs(65535 - salida_ip.get(x,y)));
                            }
                    }
                }
            imp.updateAndDraw();
            IJ.wait(500);           
            } 

            else if (resolution != 16){
                IJ.run("Invert", "stack");
                IJ.wait(500);
            }
        }

        ImagePlus imp = IJ.getImage(); 
        ImageStack stack = imp.getStack(); 
        Roi cuadradointeres;                            // variable tipo Roi
        Rectangle roi;                              // Region of interest multifractal
        ImageProcessor salida_ip = stack.getProcessor(1);
        Rectangle fullimage = salida_ip.getRoi();

        if (resolution == 16){
            //Rectangle r = salida_ip.getRoi();
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
            IJ.run("Invert", "stack");
            IJ.wait(500);
        }

        long time_start, time_end;
        time_start = System.currentTimeMillis();


    	/** ---------Start of computational analysis---------------------------------- */         
        int epsilonmax  =  (int)(floor(log(w)/log(2))); 
        int pixelintro;
        int pixelused;
        GenericDialog gd2 = new GenericDialog("Select min size of box");
        gd2.addNumericField("Minimum size of box (in pixels: ", 1, 0);
        gd2.showDialog();
            if (gd2.wasCanceled()) return;
        pixelintro       = (int)gd2.getNextNumber();
        
        pixelused = (int)pow(2,(int)(log(pixelintro)/log(2)));
        IJ.log("The minimum size of box that will be used is: "+ pixelused);
        
        double sum = 0;
        for (int z=1; z<stack.getSize()+1; z++){ 
            ImageProcessor salidaux_ip = stack.getProcessor(z);
            sum += salidaux_ip.getStatistics().area*salidaux_ip.getStatistics().mean;
        }

        double P = 0;

        epsilonmax = (int)((log(fullimage.width/pixelused))/log(2));
        double[] epsilon     = new double[epsilonmax+1]; 
        double[][] Xnum      = new double[epsilonmax+1][q.length];
        double[][] denomMu       = new double[epsilonmax+1][q.length];
        double[][] numalpha      = new double[epsilonmax+1][q.length];
        double[][] numf          = new double[epsilonmax+1][q.length];

    	/** ---------Box Counting Loop---------------------------------------------------- */         
        for (int n = 0; n<epsilon.length; n++){
            IJ.log(""+(epsilon.length-n));
            epsilon[n]         = 1/pow(2,n);
            int numpixel       = fullimage.width/(int)pow(2,n);
            int box            = 0;
            double boxnumbers  = pow(2,n*2)*pow(2,n);
            double[] Parchive  = new double[(int)boxnumbers];


            for (int z=1; z<stack.getSize()+1; z=z+numpixel){
                for (int y=fullimage.y; y<fullimage.y+fullimage.height; y=y+numpixel){
                    for (int x=fullimage.x; x<fullimage.x+fullimage.width; x=x+numpixel){
                        P = 0; 
                        for (int zz=z; zz<z+numpixel; zz++){
                            salida_ip = stack.getProcessor(zz); 
                            cuadradointeres = new Roi(x, y, numpixel, numpixel);
                            salida_ip.setRoi(cuadradointeres);
                            roi = salida_ip.getRoi();

                            for (int yy=roi.y; yy<(roi.y+roi.height); yy++){
                                for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
                                    P += salida_ip.get(xx,yy)/sum;
                                }
                            }
                        
                        }
                        Parchive[box] = P;
                        box++;

                    	/** ---------Xnum---------------------------------- */         
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
            /** ---------Mu---------------------------------- */         
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
            //Rectangle r = salida_ip.getRoi();
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
            IJ.run("Invert", "stack");
            IJ.wait(500);
        }

    	/** ---------Handle Results---------------------------------- */         
        ToolMultifr result = new ToolMultifr();
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
