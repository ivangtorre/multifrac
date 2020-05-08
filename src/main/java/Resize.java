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
import static java.lang.Math.floor;
import static java.lang.Math.log;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.process.Blitter;
import ij.process.ImageProcessor;

/**
 * Handles resize or cut image if necessary
 * 
 * @author Ivan G Torre
 */

public class Resize{
	ImageProcessor ip;
	
	public ImagePlus checkandresize(){
		/** ----------Check size ----------------------------- */
		double woriginal = ip.getWidth();
		double w = floor(log(woriginal)/log(2));
		double horiginal  = ip.getHeight();	
		double h = floor(log(horiginal)/log(2));
		w = h = min(pow(2,h),pow(2,w));	
		int choose;
		
		/** ----------Dialog ----------------------------- */
		ImagePlus salida = NewImage.createByteImage("Analyzed image", (int)w, (int)h, 1, NewImage.FILL_BLACK); 
		if (woriginal!=w || horiginal!=h){
			GenericDialog gd = new GenericDialog("Cut or resize image");
			String[] titleArray2 = new String[2];
			gd.addMessage("The image has not the properly size, choose what to do");
			titleArray2[0] = "Cut image";
			titleArray2[1] = "Resize image";
			gd.addChoice("  ", titleArray2, titleArray2[0]);
			gd.showDialog();
					if (gd.wasCanceled()) return salida;
			choose = gd.getNextChoiceIndex();
			
			/** ----------Cut ----------------------------- */
			if (choose == 0){
				Roi subpart = new Roi(0, 0, w, h);
				ip.setRoi(subpart);
				IJ.run("Copy", "");
				salida = NewImage.createByteImage("Analyzed image", (int)w, (int)h, 1, NewImage.FILL_BLACK);
				salida.show();			
				IJ.run("Paste", "");	
				salida.show();
			}
			
			/** ----------Resize ----------------------------- */
			else if (choose == 1){	
				salida = NewImage.createByteImage("Analyzed image", (int)woriginal, (int)horiginal, 1, NewImage.FILL_BLACK);
				ImageProcessor salida_aux = salida.getProcessor(); 				
				salida_aux.copyBits(ip, 0, 0, Blitter.COPY);					
				salida.show();				
				salida.updateAndDraw();	
				IJ.run("Size...", "width="+w+" height="+h+" average interpolation=Bilinear");
			}
		}
		
		/** ----------Do nothing ----------------------------- */
		else if (woriginal==w && horiginal==h) {
			salida = NewImage.createByteImage("Analyzed image", (int)w, (int)h, 1, NewImage.FILL_BLACK);
			ImageProcessor salida_aux = salida.getProcessor(); 	
			salida_aux.copyBits(ip,0,0,Blitter.COPY);			
			salida.show();					
			salida.updateAndDraw();
		}
		return salida;
	}
}
