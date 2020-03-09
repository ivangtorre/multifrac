// Ivan Gonzalez Torre
// Libraries import
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
import static java.lang.Math.*;
import java.util.ArrayList;
import java.lang.*;

public class _3D_Structure_Function implements PlugInFilter {

protected ImageStack stackaux;

	public int setup(String arg, ImagePlus impaux) {
		stackaux = impaux.getStack();
		if (impaux.getProcessor().isInvertedLut()){
			IJ.run("Invert LUT");
		}
		return DOES_8G+DOES_16+DOES_32+STACK_REQUIRED;				// Suports gray scale image
	}

	public void run(ImageProcessor ip) {						// Process image and save in IP

		//DIALOG Tau
		taudialog tauconstruct = new taudialog();
		int[] lista_tau;
		lista_tau = tauconstruct.gettau();

        //INPUT DIALOG Q VALUES ********************************************************************************
         qdialog ask_q = new qdialog();
         double[] q = ask_q.getq();
	

		////////////////////////////////////////////////
		////////////// DUPLICATE IMAGE /////////////////
		////////////////////////////////////////////////
		IJ.run("Options...", "iterations=1 black count=1"); 
		IJ.run("Duplicate...","title=Analyzed_image.tif duplicate");		// Duplicate image
		IJ.run("Options...", "iterations=1 black count=1");			// Check not invertign LUT

		ImagePlus imp = IJ.getImage(); 
		ImageStack stack = imp.getStack(); 
		Roi cuadradointeres;							// variable tipo Roi
		Rectangle roi;								// Region of interest multifractal
		ImageProcessor salida_ip = stack.getProcessor(1);
		Rectangle fullimage = salida_ip.getRoi();

		// TIME //
		long time_start, time_end;
		time_start = System.currentTimeMillis();


		// GETTING THE MATRIX FROM IMAGE
		double infinito = pow(0,-1);
		int tau, m = fullimage.height, n = fullimage.width;			// Matrix is (m rows x n columns) size
		double fxy_tau, autocorr_xy;		
		
		double[][][] matriz = new double[n][m][stack.getSize()];
		for (int z=1; z<stack.getSize()+1; z++){
			ImageProcessor salidaux_ip = stack.getProcessor(z);						
			for (int y=fullimage.y; y<(fullimage.y+m); y++){			
				for (int x=fullimage.x; x<(fullimage.x+n); x++){
					matriz[x][y][z-1] = salidaux_ip.get(x,y);
				}
			}
		}
		
		////////////////////////////////////////////////////
		//////////////MAIN OPERATIONS///////////////////////
		////////////////////////////////////////////////////
		IJ.log("Calculating...");
		double[][] Mq = new double[lista_tau.length][q.length];		
		for (int qposition = 0; qposition < q.length; qposition++){
			IJ.log("" + (q.length - qposition));
			for (int tauposition = 0; tauposition < lista_tau.length; tauposition++){
				tau = lista_tau[tauposition];
				// For a given q and tau do autocorrelation
				ArrayList<Double> lista = new ArrayList<Double>();
				for (int z=0; z<stack.getSize(); z++){
					for (int y = 0; y<m; y++){
						for (int x = 0; x<n; x++){
							fxy_tau = matriz[(x+tau)%n][y][z] + matriz[(x-tau+n)%n][y][z];
							fxy_tau += matriz[x][(y+tau)%m][z] + matriz[x][(y-tau+m)%m][z];
							fxy_tau += matriz[x][y][(z+tau)%stack.getSize()] + matriz[x][y][(z-tau+stack.getSize())%stack.getSize()];
							fxy_tau = fxy_tau/6;
							autocorr_xy = pow(abs(fxy_tau - matriz[x][y][z]), q[qposition]);
							if (autocorr_xy == infinito) {autocorr_xy = 0;}
							lista.add(autocorr_xy);
						}
					}
				}
				double media, sum = 0.0;
				for (int i = 0; i < lista.size(); i++) {
					sum += lista.get(i);
				}
				Mq[tauposition][qposition] = sum/lista.size();
			}
		}
		

		////////////////////////////////////////////////////
		///////// FIGURES AND RESULTS///////////////////////
		////////////////////////////////////////////////////
		resultestructure result = new resultestructure();
		result.Mq = Mq;
		result.lista_tau = lista_tau;
		result.lista_q = q;
		result.out();

		// TOTAL TIME
		time_end = System.currentTimeMillis();
		IJ.log("");
		IJ.log("The task has taken "+ ( time_end - time_start )/1000 +" seconds");

	}
}

