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

public class _2D_Structure_Function implements PlugInFilter {
public ImagePlus salida;
    public int setup(String arg, ImagePlus imp) {
        if (imp.getProcessor().isInvertedLut()){
            IJ.run("Invert LUT");
        }
        return DOES_8G+DOES_16+DOES_32;                     // Suports gray scale image
        
    }

    public void run(ImageProcessor ip) {                        // Process image and save in IP

        //DIALOG Tau
        taudialog tauconstruct = new taudialog();
        int[] lista_tau;
        lista_tau = tauconstruct.gettau();

        //DIALOG Q
        //INPUT DIALOG Q VALUES ********************************************************************************
         qdialog ask_q = new qdialog();
         double[] q = ask_q.getq();
        //qdialogtau qconstruct = new qdialogtau();
        //double[] lista_q;
        //lista_q = qconstruct.getq();
    

        ////////////////////////////////////////////////
        ////////////// DUPLICATE IMAGE /////////////////
        ////////////////////////////////////////////////
        int resolution = ip.getBitDepth();
        ImagePlus salida = NewImage.createImage("Analyzed Image", (int)ip.getWidth(), (int)ip.getHeight(), 1,resolution, NewImage.FILL_BLACK);
        ImageProcessor salida_aux = salida.getProcessor();      // Processing
        salida_aux.copyBits(ip,0,0,Blitter.COPY);           // Copy pixels  
        salida.show();                          // Show
        salida.updateAndDraw();     
        ImageProcessor salida_ip = salida.getProcessor();
        salida_ip = salida.getProcessor();              // Processing the image after the changes
        Rectangle fullimage = salida_ip.getRoi();           // ROI full image
    
        // TIME //
        long time_start, time_end;
        time_start = System.currentTimeMillis();


        // GETTING THE MATRIX FROM IMAGE
        double infinito = pow(0,-1);
        int tau, m = fullimage.height, n = fullimage.width;         // Matrix is (m rows x n columns) size
        double fxy_tau, autocorr_xy;     
        
        double[][] matriz = new double[n][m];                       
        for (int y=fullimage.y; y<(fullimage.y+m); y++){            
            for (int x=fullimage.x; x<(fullimage.x+n); x++){
                matriz[x][y] = salida_ip.get(x,y);
            }
        }
        
        ////////////////////////////////////////////////////
        //////////////MAIN OPERATIONS///////////////////////
        ////////////////////////////////////////////////////
        IJ.log("Calculating...");
        double[][] Mq = new double[lista_tau.length][q.length];
        for (int qposition = 0; qposition < q.length; qposition++){
            IJ.log(""+(q.length - qposition)); 
            //q = q[qposition];
            for (int tauposition = 0; tauposition < lista_tau.length; tauposition++){
                tau = lista_tau[tauposition];
                // For a given q and tau do autocorrelation
                ArrayList<Double> lista = new ArrayList<Double>();
                for (int y = 0; y<m; y++){
                    for (int x = 0; x<n; x++){
                        fxy_tau = matriz[(x+tau)%n][y]+matriz[(x-tau+n)%n][y]+matriz[x][(y+tau)%m]+matriz[x][(y-tau+m)%m];
                        fxy_tau = fxy_tau/4;
                        autocorr_xy = pow(abs(fxy_tau - matriz[x][y]),q[qposition]);
                        if (autocorr_xy == infinito) {autocorr_xy = 0;}
                        lista.add(autocorr_xy);
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

