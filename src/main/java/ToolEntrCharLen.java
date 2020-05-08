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
 * Tool that computes Configuration entropy and characteristic length u
 * using Gliding Box sliding method on 2D Gray images
 * 
 * @author Ivan G Torre
 */

public class ToolEntrCharLen{
	ImageProcessor ip;
	String str;
	String path;

	public void out() {

		int woriginal  = ip.getWidth();
		int horiginal  = ip.getHeight();
		int lmax = min(woriginal,horiginal);
		double[] lsize = new double[(int)lmax];
		//double[] lacunarity = new double[(int)lmax];
		double[] entropy = new double[(int)lmax+1];
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
			double numblackpixel;
			double total_pixels;					
			int total_pixel = n*n;
			
			for (int y=fullimage.y; y<fullimage.y+fullimage.height-numpixel+1; y++){
				for (int x=fullimage.x; x<fullimage.x+fullimage.width-numpixel+1; x++){	
					cuadradointeres = new Roi(x, y, numpixel, numpixel);
					ip.setRoi(cuadradointeres);
					roi = ip.getRoi();
					numblackpixel=0;
					total_pixels = 0;

					for (int yy=roi.y; yy<(roi.y+roi.height); yy++){		
						for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
							total_pixels += 1;
							if (ip.get(xx,yy) == 0){numblackpixel += +1;}
						}
					}
					measure_box[counter] = total_pixels - numblackpixel;
					counter += 1;
				}
			}

	    	/** ---------Probability mass vector--------------------------------- */         
			double[] freq_distribution = new double[(int)total_pixel + 1];
			for (int pos = 0; pos<measure_box.length; pos++){
				freq_distribution[(int)measure_box[pos]] += 1;
			}

			for (int dum = 0; dum<freq_distribution.length; dum++){
				freq_distribution[dum] = (double)freq_distribution[dum]/(double)boxnumbers;
			}
			
	    	/** ---------Compute Configuration Entropy--------------------------------- */         
			double local_entropy = 0;	 
			//double z1 = 0;
			//double z2 = 0;		
			for (int m = 0; m<freq_distribution.length; m++){
				if (freq_distribution[m] > 0){
					//z1 += pow(m,1) * (double)freq_distribution[m];
					//z2 += pow(m,2) * (double)freq_distribution[m];
					local_entropy += -(double)freq_distribution[m] * (double)log(freq_distribution[m]);
				}
			}
	    	/** ---------Normalize entropy--------------------------------------------------- */         
			local_entropy = (double)local_entropy/(double)log(n*n + 1);
			entropy[lsize_iter] = local_entropy;				// Entropy vector	
			//lacunarity[lsize_iter] = (double)z2/(pow(z1,2)); 		// Lacunarity vector	
			lsize[lsize_iter] = n;
			lsize_iter += 1;
		}
		
    	/** ---------Get Characteristic length--------------------------------------------- */         
		double max_entropy = 0;
		double char_length = 0;
		for (int dumm = 1; dumm < entropy.length; dumm++){
			if (entropy[dumm] > max_entropy){
				max_entropy = entropy[dumm];
				char_length = lsize[dumm];
			}
		}
		IJ.log("char_length " + char_length);
		

    	/** ---------Handle Results------------------------------------------------------ */         
		String str = new String(" 2D Multifractal Gliding method");

    	/** ---------Plot Configuration entropy---------------------------------------- */
		PlotWindow.noGridLines = false;
		Plot plot_entropy = new Plot("H*(l)", "l", "H*(l)");
		plot_entropy.add("line", lsize, entropy);
		plot_entropy.setLimits(1, lsize[lsize.length-1], 0, 1.3 );
		plot_entropy.setLineWidth(2);
		plot_entropy.show();

    	/** ---------Save Configuration entropy plot---------------------------------- */         
		ImagePlus imp_plot1 = plot_entropy.getImagePlus();
		IJ.saveAs(imp_plot1, "tif", path + "Conf_entropy"+ str);
		//imp_plot1.close();
		

    	/** ---------Results table Configuration Entropy---------------------------------- */         
 		ResultsTable table = new ResultsTable();
 		for(int dumm_c = 0; dumm_c < lsize.length; dumm_c++)    {
			table.incrementCounter();
			table.addValue("lsize", IJ.d2s(lsize[dumm_c], 2));
			table.addValue("conf_entropy",IJ.d2s(entropy[dumm_c], 8));
		}

		table.setPrecision(8);
		table.updateResults();
		table.setPrecision(8);

		table.show("Results");
		IJ.selectWindow("Results");
		IJ.save(path+"conf_entrop.txt");
		//IJ.run("Close");
		
 		ResultsTable table2 = new ResultsTable();
		table2.incrementCounter();
		table2.addValue("Characteristic length", IJ.d2s(char_length, 2));
		//table2.addValue("Lacunarity", IJ.d2s(slope, 5));
		//table2.addValue("Lacunarity error", IJ.d2s(devslope, 5));

		table2.setPrecision(8);
		table2.updateResults();
		table2.setPrecision(8);

		table2.show("Results2");
		IJ.selectWindow("Results2");
		IJ.save(path+"char_length.txt");
		//IJ.run("Close");
		//IJ.selectWindow("Results");
		//IJ.run("Close");
		

		// LACUNARITY PLOT
		/*PlotWindow.noGridLines = false; 							// draw grid lines
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
		//cf       = fit1.getcf();							// contructor	
		//plot(cf);									// plot

		// Print slope and its error
		//IJ.log("");
		//IJ.log("Box Counting Lacunarity:");
		//IJ.log("");
		//IJ.log("D = "+IJ.d2s(slope,5)+" +/- "+IJ.d2s(devslope,5));

*/
		//IJ.log(str);
		//table1.setPrecision(8);
		//table1.updateResults();
		//table1.setPrecision(8);
		//table1.show("Results");
		//IJ.selectWindow("Results");
		//IJ.save(path+"Results.txt");
	}
}
