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

import java.awt.*;
import ij.*;
import ij.gui.*;
import ij.process.*;
import ij.measure.*;
import static java.lang.Math.*;

/**
 * Tool that computes Lacunarity analysis using Gliding Box sliding method
 * on 2D Gray images
 * 
 * @author Ivan G Torre
 */

public class ToolLacun{
	ImageProcessor ip;
	String str;
	String path;

	public void out() {

		int woriginal  = ip.getWidth();					
		int horiginal  = ip.getHeight();
		int lmax = min(woriginal,horiginal);
		double[] lsize = new double[(int)lmax];
		double[] lacunarity = new double[(int)lmax];
		Roi cuadradointeres;
		Rectangle roi;				
		Rectangle fullimage = ip.getRoi();

    	/** ---------Gliding Box Loop---------------------------------- */         
		int lsize_iter = 0;
		for (int n = 1; n<=lmax; n++){
			IJ.log(""+(lmax-n));
			int counter = 0;
			int numpixel = (n);
			//int box      = 0;
			double boxnumbers  = (woriginal-n+1)*(horiginal-n+1);
			double[] measure_box  = new double[(int)boxnumbers];
			//double[] measure_box_counting  = new double[(int)boxnumbers];
			//double numblackpixel;
			double total_pixels=0;					

			for (int y=fullimage.y; y<fullimage.y+fullimage.height-numpixel+1; y++){
				for (int x=fullimage.x; x<fullimage.x+fullimage.width-numpixel+1; x++){	
					cuadradointeres = new Roi(x, y, numpixel, numpixel);
					ip.setRoi(cuadradointeres);
					roi = ip.getRoi();
					int mass=0;
					total_pixels = 0;
					
			    	/** ---------Count pixels along roi---------------------------------- */         
					for (int yy=roi.y; yy<(roi.y+roi.height); yy++){		
						for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
							total_pixels += 255;
							mass += ip.get(xx,yy);
						}
					}
					
					measure_box[counter] = mass;
					counter += 1;

				}
			}

	    	/** ---------Probability mass vector--------------------------------- */         
			double[] freq_distribution = new double[(int)total_pixels + 1];	
			for (int pos = 0; pos<measure_box.length; pos++){
				freq_distribution[(int)measure_box[pos]] += 1;
			}

			for (int dum = 0; dum<freq_distribution.length; dum++){
				freq_distribution[dum] = (double)freq_distribution[dum]/(double)boxnumbers;
			}
			
	    	/** ---------Compute Lacunarity--------------------------------- */         
			double z1 = 0;
			double z2 = 0;		
			for (int m = 0; m<freq_distribution.length; m++){
				if (freq_distribution[m] > 0){
					z1 += pow(m,1) * (double)freq_distribution[m];
					z2 += pow(m,2) * (double)freq_distribution[m];
				}
			}
			lacunarity[lsize_iter] = (double)z2/(pow(z1,2));
			lsize[lsize_iter] = n;
			lsize_iter += 1;	
		}
		
    	/** ---------Handle Results---------------------------------- */         
		String str = new String(" 2D Multifractal Gliding method");
		
    	/** ---------Plot Lacunarity---------------------------------- */         
		PlotWindow.noGridLines = false;
		Plot plot_lacunarity = new Plot("lacunarity(l)", "l", "Lacunarity");
		plot_lacunarity.add("line", lsize, lacunarity);
		plot_lacunarity.setLineWidth(2);
		plot_lacunarity.show();

    	/** ---------Save Lacunarity plot---------------------------------- */         
		ImagePlus imp_plot0 = plot_lacunarity.getImagePlus();
		IJ.saveAs(imp_plot0, "tif", path + "Lacunarity"+ str);
		//imp_plot0.close();


    	/** ---------Fit Lacunarity---------------------------------- */         
		double[] x = new double[lmax];
		double[] y = new double[lmax];
		for(int counter = 0; counter < lmax; counter++){
			IJ.log("" + lsize[counter]);
            x[counter] = log(lsize[counter]);
            y[counter] = log(lacunarity[counter]);
        	}

		Linefit fit1 = new Linefit();
		double slope, devslope;
		//CurveFitter cf;
		fit1.x   = x;
		fit1.y   = y;
		slope    = fit1.getslope();
		devslope = fit1.getdevslope();
		

    	/** ---------Results table lacunarity---------------------------------- */         
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
		//IJ.run("Close");

		
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
		//IJ.run("Close");
		//IJ.selectWindow("Results");
		//IJ.run("Close");

	}
}
