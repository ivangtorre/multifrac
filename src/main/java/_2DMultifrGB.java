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
import ij.measure.*;
import java.io.File;
import static java.lang.Math.*;

/**
 * Tool that computes multifractal analysis using Gliding Box sliding method
 * on 2D Gray images
 * 
 * @author Ivan G Torre
 */

public class _2DMultifrGB implements PlugInFilter {
    public String base_path;
    public String imagename;    
    public int setup(String arg, ImagePlus imp) {
        //if (imp.getProcessor().isInvertedLut()){
        //  IJ.run("Invert LUT");
        //}
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
				"	Torre, I.G., Heck R.J. & Tarquis, A.M. (2020).\n" +
				"   MULTIFRAC: An ImageJ plugin for multiscale characterization\n" + 
				"   of 2D and 3D stack images.");
    	
    	/** ----------Dialogs and options------------- */      
    	Qdialog ask_q = new Qdialog();
    	double[] q = ask_q.getq();

    	/** ---------Check size of image------------- */         
        double woriginal  = ip.getWidth(); 
        double horiginal  = ip.getHeight();
        double maxsize    = min(woriginal,horiginal);

        GenericDialog gd = new GenericDialog("Select min and max size of box (in pixels)");
        gd.addNumericField("minsize: ", 1, 0);
        gd.addNumericField("maxsize: ", maxsize, 0);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        int nmin       = (int)gd.getNextNumber();
        int nmax       = (int)gd.getNextNumber();

        int resolution = ip.getBitDepth();
        ImagePlus salida = NewImage.createImage("Analyzed Image", (int)woriginal, (int)horiginal, 1, resolution, NewImage.FILL_BLACK);
        ImageProcessor salida_aux = salida.getProcessor();
        salida_aux.copyBits(ip,0,0,Blitter.COPY);
        salida.show();
        salida.updateAndDraw();     


    	/** ---------Choose which color is max value------------- */         
        ImageProcessor salida_ip = salida.getProcessor();
        GenericDialog gd1   = new GenericDialog("Process image");           // Create dialog
        String[] titleArray2    = new String[2];
        gd1.addMessage("Maximum measure will be for the white. Do you want to analyze this image or the inverted one?");
        titleArray2[0] = "This image";
        titleArray2[1] = "Inverted image";
        gd1.addChoice("  ", titleArray2, titleArray2[0]);
        gd1.showDialog();
        if (gd1.wasCanceled()) return;
        int choose;
        choose = gd1.getNextChoiceIndex();
        if (choose == 0){
        }
        else if (choose == 1){// Then invert image
            salida_ip = salida.getProcessor();      // Processing the image after the changes
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
                salida.updateAndDraw();
                IJ.wait(500);
            }
        }

    	/** ---------Handle output path---------------------------------- */         
        String save_path_pred;
        String[] parts = imagename.split("\\.");
        String part1 = parts[0]; // 004
        base_path = base_path + part1 + "/";
        if (choose == 0){save_path_pred = base_path + "gliding/white_max_value/min_" + nmin + "px/";} 
        else {save_path_pred = base_path + "gliding/black_max_value/min_" + nmin + "px/";}
        String path;
        path = save_path_pred;

    	/** ---------Check if path exists and if not, create it----------------------------- */         
        File wdir = new File(path);
        if (!wdir.isDirectory() || !wdir.exists() || !wdir.canRead()) {
            wdir.mkdirs();
            IJ.log("Path created.");}


        salida_ip = salida.getProcessor();
        Rectangle fullimage = salida_ip.getRoi();
        if (resolution == 16){
            Rectangle r = salida_ip.getRoi();
                for (int y=r.y; y<(r.y+r.height); y++){
                    for (int x=r.x; x<(r.x+r.width); x++){
                        salida_ip.set(x, y, (int)abs(65535 - salida_ip.get(x,y)));
                    }
            }
        salida.updateAndDraw();
        IJ.wait(500);
        } 

        else if (resolution != 16){}
        
        long time_start, time_end;
        time_start = System.currentTimeMillis();

    	/** ---------Start of computational analysis---------------------------------- */         
        double sum = salida_ip.getStatistics().area*salida_ip.getStatistics().mean;
        double[][] matriz = new double[fullimage.width][fullimage.height];
        for (int y=fullimage.y; y<(fullimage.y+fullimage.height); y++){
            for (int x=fullimage.x; x<(fullimage.x+fullimage.width); x++){
                matriz[x][y] = salida_ip.get(x,y)/sum;
            }
        }

        double[] epsilon    = new double[(int)maxsize+1];
        double[][][] Xnum       = new double[(int)maxsize+1][q.length][3];         
        double deltavalue   = 0.01;     
        double [] delta     = {0,-deltavalue,deltavalue};
        double P = 0;
        Roi cuadradointeres;
        Rectangle roi;

    	/** ---------Gliding Box sliding method Loop---------------------------------- */         
        for (int n = nmin; n<=nmax; n++){
            IJ.log(""+(nmax-n));                            // Show bucles left 
            epsilon[n]   = (double)(n)/(double)maxsize;
            int numpixel = (n);                        
            int box      = 0;
            double boxnumbers  = (woriginal-n+1)*(horiginal-n+1);
            double[] Parchive  = new double[(int)boxnumbers];

            // Creating ROIs along image
            for (int y=fullimage.y; y<fullimage.y+fullimage.height-numpixel+1; y++){
                for (int x=fullimage.x; x<fullimage.x+fullimage.width-numpixel+1; x++){
                    cuadradointeres = new Roi(x, y, numpixel, numpixel);    // Creating size of ROI
                    salida_ip.setRoi(cuadradointeres);              // Fix ROI in the image
                    roi = salida_ip.getRoi();               // Processing ROI
                    P=0;
                
                    // Counting pixels along roi 
                    for (int yy=roi.y; yy<(roi.y+roi.height); yy++){        
                        for (int xx=roi.x; xx<(roi.x+roi.width); xx++){
                            P += matriz[xx][yy];
                        }
                    }
                    if (P>0) {
                        Parchive[box] = P;
                        box++;
                    	/** ---------Xnum---------------------------------- */         
                        for (int qposition = 0; qposition<q.length; qposition++){
                            for (int deltaposition = 0; deltaposition<3; deltaposition++){
                                Xnum[n][qposition][deltaposition] += pow(P,q[qposition]+delta[deltaposition]);
                            }
                        }   
                    }
                }
            }
            for (int qposition = 0; qposition<q.length; qposition++){
                for (int deltaposition = 0; deltaposition<3; deltaposition++){
                    Xnum[n][qposition][deltaposition] = Xnum[n][qposition][deltaposition]/box;
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
            IJ.run(salida, "Invert", "");                       // Invert image. Count black mass
            salida.updateAndDraw();
            IJ.wait(500);
        }

        
    	/** ---------Handle Results---------------------------------- */         
        double minx0 = 0;
        double maxx0 = 0;
        double[] x = new double[nmax-nmin+1];
        for(int n = 0; n <= nmax-nmin; n++){
            x[n] = log(epsilon[n+nmin]);
            if (x[n]>maxx0){maxx0=x[n];}
            if (x[n]<minx0){minx0=x[n];}
        }

        double miny0 = 0;
        double maxy0 = 0;   
        double[][][] y = new double[nmax-nmin+1][q.length][3]; 
        for (int n = 0; n<=nmax-nmin; n++){
            for (int qposition = 0; qposition<q.length; qposition++){
                for (int deltaposition = 0; deltaposition<3; deltaposition++){
                    y[n][qposition][deltaposition] = log(Xnum[n+nmin][qposition][deltaposition]);
                    if (y[n][qposition][deltaposition]>maxy0){maxy0=y[n][qposition][deltaposition];}
                    else if (y[n][qposition][deltaposition]<miny0){miny0=y[n][qposition][deltaposition];}
                }       
            }
        }       
        
        // CURVE FITTING (type of fitType: STRAIGHT_LINE=0, POLY2=1, POLY3=2, POLY4=3, EXPONENTIAL=4, POWER=5,...)  
        double[] TauQ       = new double[q.length]; 
        double[][] TauQaux      = new double[q.length][3];
        double[] desvTauQ   = new double[q.length];
        double[][] desvTauQaux  = new double[q.length][3];
        double[] relerror   = new double[q.length];
        double[] Dq         = new double[q.length];
        double[] desvDq     = new double[q.length];
        double[] alphaq     = new double[q.length];
        double[] desvalphaq     = new double[q.length];
        double[] fq         = new double[q.length];
        double[] desvfq     = new double[q.length];
        double[] ybar       = new double[y.length];

        // TauQ and eror
        for (int qposition = 0; qposition<q.length; qposition++){
            double[] yaux = new double[y.length];
            for (int deltaposition = 0; deltaposition<3; deltaposition++){
                for(int n = 0; n < y.length; n++){
                    yaux[n] = y[n][qposition][deltaposition];
                }
            
                Linefit fit1 = new Linefit();
                fit1.x   = x;
                fit1.y   = yaux;
                TauQaux[qposition][deltaposition]       = fit1.getslope() - 2;
                desvTauQaux[qposition][deltaposition]   = fit1.getdevslope();
                if (deltaposition == 0){
                    TauQ[qposition]         = fit1.getslope() - 2;
                    desvTauQ[qposition] = fit1.getdevslope();
                }
            }
            relerror[qposition] = desvTauQ[qposition]/abs(TauQ[qposition]);
            if (q[qposition] == 1){
                relerror[qposition] = (relerror[qposition-1]+relerror[qposition+1])/2;
            }
        }

        // Dq and error
        for (int qposition = 0; qposition<q.length; qposition++){
            if (q[qposition] != 1){
                Dq[qposition]     = TauQ[qposition]/(q[qposition] - 1);
                desvDq[qposition] = desvTauQ[qposition]/abs(q[qposition] - 1);
            }   
        }
        for (int qposition = 0; qposition<q.length; qposition++){
            if (q[qposition] == 1){
                Dq[qposition]     = (Dq[qposition-1]+Dq[qposition+1])/2;
                desvDq[qposition] = (desvDq[qposition-1]+desvDq[qposition+1])/2;
            }
        }

        
        /** --------PLOTS------------------------------------------------------ */ 
        String str = new String(" 2D Multifractal Gliding method");

		/** --------Plot Xq-epsilon------------------------------------------------------ */ 
        Plot plot0 = new Plot("X(q,epsilon)"+str, "ln(epsilon)", "ln(X(q))");
        PlotWindow.noGridLines = false;
        plot0.setLimits(minx0-0.5, maxx0 + 0.5,miny0-1,maxy0+1);
        for(int qposition = 0; qposition < q.length; qposition++){
            for(int n = 0; n < y.length; n++){
                ybar[n] = y[n][qposition][0];
            }
        plot0.addPoints(x, ybar, Plot.LINE);
        }
        plot0.show();
        ImagePlus imp_plot0 = plot0.getImagePlus();
        IJ.saveAs(imp_plot0, "Jpeg", path + "X(q,epsilon)"+ str);
        imp_plot0.close();
        IJ.selectWindow("X(q,epsilon)"+ str);
        
        
		/** --------Plot Tauq------------------------------------------------------ */ 		
        PlotWindow.noGridLines = false;
        Plot plot1 = new Plot("Tau(q)"+str, "q", "TauQ");
		plot1.add("line", q, TauQ);
        plot1.setLimits((q[0])-0.5, q[q.length -1]+0.5, TauQ[0]-desvTauQ[q.length -1]-0.5, TauQ[q.length -1]+desvTauQ[0]+0.5);
        plot1.setLineWidth(2);
        plot1.addErrorBars(desvTauQ);
        plot1.setColor(Color.red);
        plot1.show();
        ImagePlus imp_plot1 = plot1.getImagePlus();
        IJ.saveAs(imp_plot1, "Jpeg", path + "TauQ vs q"+ str);
        imp_plot1.close();

        
		/** --------Plot Dq------------------------------------------------------ */
        PlotWindow.noGridLines = false;
        Plot plot2 = new Plot("D(q)"+str, "q", "Dq");
		plot2.add("line",q, Dq);
        plot2.setLimits((q[0]) -0.5, q[q.length -1] + 0.5, Dq[q.length -1] - desvDq[q.length -1] - 0.5, Dq[0] + desvDq[0] + 0.5);
        plot2.setLineWidth(2);
        plot2.addErrorBars(desvDq);
        plot2.setColor(Color.red);
        plot2.show();
        ImagePlus imp_plot2 = plot2.getImagePlus();
        IJ.saveAs(imp_plot2, "Jpeg", path + "Dq vs q"+ str);
        imp_plot2.close();


        // We estimate alpha(q) with central differences methods. 
        // alpha(q) = dTau(q)/dq;
        double maxx         = 0;
        double minx         = 10;
        double maxdesvalphaq    = 0;
        double alpha1       = 0;
        double err_alpha1   = 0.0;
        double alphamin     = 100;
        double err_alphamin     = 0;
        int qposition_alphamin  = 0;

        double alphamax     = -100;
        double err_alphamax     = 0;
        int qposition_alphamax  = 0;

        for (int qposition = 0; qposition<q.length; qposition++){
            alphaq[qposition] = (TauQaux[qposition][2]-TauQaux[qposition][1])/((q[qposition]+deltavalue) - (q[qposition]-deltavalue));
            if (alphaq[qposition]<minx){minx=alphaq[qposition];}
            if (alphaq[qposition]>maxx){maxx=alphaq[qposition];}
            desvalphaq[qposition] = abs((alphaq[qposition])*(relerror[qposition]));
            if (desvalphaq[qposition]>maxdesvalphaq){maxdesvalphaq=desvalphaq[qposition];}
            if (q[qposition] == 1){alpha1 = alphaq[qposition];
                        err_alpha1 = desvalphaq[qposition];}
            if (alphaq[qposition] < alphamin){alphamin = alphaq[qposition];
                                err_alphamin = desvalphaq[qposition];
                                qposition_alphamin = qposition;}
            if (alphaq[qposition] > alphamax){alphamax = alphaq[qposition];
                                err_alphamax = desvalphaq[qposition];
                                qposition_alphamax = qposition;}

        }

		/** --------Plot alfa(q)---------------------------------- */
        PlotWindow.noGridLines = false;                         // draw grid lines
        Plot plot3 = new Plot("alpha(q)"+str, "q", "alpha(q)");
        plot3.add("line", q, alphaq);
        plot3.setLimits(q[0] -0.5, q[q.length -1] + 0.5, alphaq[q.length -1] - 0.5 - maxdesvalphaq, alphaq[0] + 0.5 +  maxdesvalphaq );
        plot3.setLineWidth(2);                              // Line width
        plot3.addErrorBars(desvalphaq);                             // Errors bars

        
        // Calculate f(alpha)
        double maxy         = 0;
        double miny         = 10;
        double maxdesvfq    = 0;

        for (int qposition = 0; qposition<q.length; qposition++){
            fq[qposition] = (alphaq[qposition])*(q[qposition]) - TauQ[qposition];
            if (fq[qposition] < miny){miny = fq[qposition];}
            if (fq[qposition] > maxy){maxy = fq[qposition];}
            desvfq[qposition] = (desvalphaq[qposition])*(abs(q[qposition])) + desvTauQ[qposition];
            if (desvfq[qposition]>maxdesvfq){maxdesvfq=desvfq[qposition];}
        }
        
		/** --------Plot f(alfa)---------------------------------- */ 		
        PlotWindow.noGridLines = false;
        Plot plot4 = new Plot("f(alpha) vs alpha"+str, "alpha", "f(alpha)");
		plot4.add("line", alphaq, fq);
        plot4.setLimits(minx - 0.5 - maxdesvalphaq,maxx + 0.5 + maxdesvalphaq,miny -0.5 - maxdesvfq,maxy + 0.5 + maxdesvfq);
        plot4.setLineWidth(2);                              // Line width
        plot4.setColor(Color.red);                          // line color
        for (int qposition = 0; qposition<q.length; qposition++){
            plot4.drawLine(alphaq[qposition]-desvalphaq[qposition],fq[qposition],alphaq[qposition]+desvalphaq[qposition],fq[qposition]);
        }

        plot4.show();
        ImagePlus imp_plot3 = plot4.getImagePlus();
        IJ.saveAs(imp_plot3, "Jpeg", path + "f(alpha) vs alpha"+ str);
        imp_plot3.close();
        
        // RESULTS TABLE
        ResultsTable table = new ResultsTable();
        for(int counter = 0; counter < q.length; counter++)    {
            table.incrementCounter();
            table.addValue("q", IJ.d2s(q[counter], 2));
            table.addValue("TauQ",IJ.d2s(TauQ[counter], 8));
            table.addValue("error TauQ",IJ.d2s(desvTauQ[counter], 8));
            table.addValue("Dq",IJ.d2s(Dq[counter], 8));
            table.addValue("error Dq",IJ.d2s(desvDq[counter], 8));
            table.addValue("alpha", IJ.d2s(alphaq[counter], 8));
            table.addValue("error alpha", IJ.d2s(desvalphaq[counter], 8));
            table.addValue("f(alpha)", IJ.d2s(fq[counter], 8));
            table.addValue("error f", IJ.d2s(desvfq[counter], 8));
        }

        for(int ep = 0; ep < x.length; ep++){
            table.incrementCounter();
            for(int qposition = 0; qposition < q.length; qposition++){
                String fs;
                fs = String.format("x %d",qposition);
                table.addValue(fs, x[ep]);
                String fs2;
                fs2 = String.format("y %d",qposition);
                table.addValue(fs2, y[ep][qposition][0]);   
            }
        }

        table.setPrecision(8);
        table.updateResults();
        table.setPrecision(8);

        table.show("Results");
        IJ.selectWindow("Results");
        IJ.save(path+"Results.txt");
        IJ.run("Close");


        // TABLE espectrum
        ResultsTable table2 = new ResultsTable();
    
        //Area_fraction
        double total_pixel = 0.0; 
        double countValue = 0.0;
        for( int d_x = 0; d_x < ip.getWidth(); d_x++ ) { 
                for( int d_y = 0; d_y < ip.getHeight(); d_y++ ) {
                total_pixel += 1; 
                        if(ip.get(d_x,d_y) == 0 ) { 
                                countValue += 1; 
                        } 
                } 
        } 


        // Write to table
        table2.incrementCounter();
        double areafraction = (total_pixel - countValue)/total_pixel;
        table2.addValue("value", IJ.d2s(areafraction, 2));
        table2.addLabel("porosity");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(alpha1, 3));
        table2.addLabel("alpha1");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(err_alpha1, 4));
        table2.addLabel("alpha1 error");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(alphamin, 5));
        table2.addLabel("alphamin");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(err_alphamin, 6));
        table2.addLabel("alphamin error");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(fq[qposition_alphamin], 7));
        table2.addLabel("f_alphamin");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(desvfq[qposition_alphamin], 8));
        table2.addLabel("f_alphamin error");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(alphamax, 9));
        table2.addLabel("alphamax");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(err_alphamax, 10));
        table2.addLabel("alphamax error");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(fq[qposition_alphamax], 11));
        table2.addLabel("f_alphamax");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(desvfq[qposition_alphamax], 12));
        table2.addLabel("f_alphamax error");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(alphamax - alphamin, 13));
        table2.addLabel("alphamax-alphamin");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(err_alphamin + err_alphamax, 14));
        table2.addLabel("alphamax-alphamin error");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(fq[qposition_alphamax] - fq[qposition_alphamin], 15));
        table2.addLabel("falphamax-falphamin");

        table2.incrementCounter();
        table2.addValue("value", IJ.d2s(desvfq[qposition_alphamax] + desvfq[qposition_alphamin], 16));
        table2.addLabel("falphamax-falphamin error");







        table2.updateResults();
        table2.show("parameteres");
        IJ.selectWindow("parameteres");
        IJ.save(path+"parameters.txt");
        IJ.run("Close");
        IJ.selectWindow("Results");
        IJ.run("Close");




        //areafraction = measures.getResult("Min");
        //table2.setPrecision(8);
        //table2.updateResults();
        //table2.setPrecision(8);

        //table2.show("Espectrum and parameters");
        //IJ.selectWindow("Espectrum and parameters");
        //IJ.save(path+"Espectrum_and_parameters.txt");



        time_end = System.currentTimeMillis();
        IJ.log("");
        IJ.log("The task has taken "+ ( time_end - time_start ) +" milliseconds");
        


    }

}
