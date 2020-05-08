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

/**
 * Tool that computes the fractal dimension of any Black and White image
 * 
 * @author Ivan G Torre
 */


public class _2DFracDimBC implements PlugInFilter {
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
		

		/** ----------Check size of image and resize if necessary----------------------------- */
		Resize checkresize = new Resize();
		checkresize.ip = ip;
		ImagePlus salida = checkresize.checkandresize();

		/** --------Convert image to B&W------------ */
		IJ.run("Options...", "iterations=1 count=1");
		IJ.run("Make Binary");

		/** --------Select which color to be counted------------ */
		ColorDialog ask_color = new ColorDialog();
        int colorchoose = ask_color.askmax();
        
		/** ---------Handle output path---------------------------------- */
		OutputPath outputpath = new OutputPath();
		outputpath.imagename = imagename;
		outputpath.base_path = base_path;
		outputpath.colorchoose = colorchoose;
		String path = outputpath.savepath();     
        
		/** ---------Save image to be analyzed---------------------------------- */
		IJ.saveAs(salida, "tif", path + "Fractal_dimension");

		/** --------Select max and min sizes for scaling------------ */		
		ImageProcessor salida_ip = salida.getProcessor();
		long time_start, time_end;
		time_start = System.currentTimeMillis();
		double w = salida_ip.getWidth();
		int epsilonmax      =  (int)(floor(log(w)/log(2)));	
		Rectangle fullimage = salida_ip.getRoi();
		Roi cuadradointeres;					
		Rectangle roi;		
		
		GenericDialog gd = new GenericDialog("Select min and max size of box (in pixels)");
		gd.addNumericField("minsize: ", 1, 0);						
		gd.addNumericField("maxsize: ", pow(2,epsilonmax), 0);	
		gd.showDialog();
		if (gd.wasCanceled()) return;		
		int pixelintromin       = (int)gd.getNextNumber();
		int pixelintromax       = (int)gd.getNextNumber();
		int pixelusedmin = (int)pow(2,(int)(log(pixelintromin)/log(2)));
		int pixelusedmax = (int)pow(2,(int)(log(pixelintromax)/log(2)));
		IJ.log("The minimum size of box that will be used is: "+ pixelusedmin);
		IJ.log("The maximum size of box that will be used is: "+ pixelusedmax);   	
		
		
		/** --------Box Counting algorithm Fractal dimension------------ */
		epsilonmax = (int)((log(fullimage.width/pixelusedmin))/log(2));
		int epsilonmin = (int)((log(fullimage.width/pixelusedmax))/log(2));
		double[] epsilon = new double[epsilonmax+1-epsilonmin];
		double[] N = new double[epsilonmax+1-epsilonmin];

		for (int n = 0; n<=epsilonmax-epsilonmin; n++){
			epsilon[n]   =  1/pow(2,n+epsilonmin);		
			int numpixel = fullimage.width/(int)pow(2,n+epsilonmin);
			
			for (int y=fullimage.y; y<fullimage.y+fullimage.height; y=y+numpixel){
				for (int x=fullimage.x; x<fullimage.x+fullimage.width; x=x+numpixel){
					cuadradointeres = new Roi(x, y, numpixel, numpixel);
					salida_ip.setRoi(cuadradointeres);
					roi = salida_ip.getRoi();
				
					/** --------Count black pixels------------ */
					bucle1:
					for (int yy=roi.y; yy<(roi.y+roi.height); yy++){		
						for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
							if (salida_ip.get(xx,yy) == 0){
								N[n]++; 
								break bucle1;
							}
						}
					}
				}
			}
		}

    	/** ---------Handle Results---------------------------------- */         
		boolean isstack = false;
        ToolFracDim result = new ToolFracDim();
		String str = new String(" 2D Monofractal Box Counting method");
		ToolFracDim.str = str;
		result.isstack = isstack;
		result.epsilon = epsilon;
		result.N = N;
		ToolFracDim.path = path;
		result.out();
		time_end = System.currentTimeMillis();
		IJ.log("The task has taken "+ ( time_end - time_start ) +" milliseconds");
	}

}

