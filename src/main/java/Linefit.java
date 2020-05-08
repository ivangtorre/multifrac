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
import ij.measure.*;
import static java.lang.Math.*;


/**
 * Compute slope and uncertainty by minimum square fit
 * 
 * @author Ivan G Torre
 */

public class Linefit {
	double[] x;
	double[] y;
	private double slope;
	private double origin;
	private double devslope;
	private static int fitType = 0;

	public CurveFitter getcf() {
		CurveFitter cf = new CurveFitter(x,y);
		cf.doFit(fitType);
		return cf;
	}

	public double getslope() {
		CurveFitter cf = new CurveFitter(x,y);
		cf.doFit(fitType);
		double[] aux = cf.getParams(); 							
		slope = aux[1];
		return slope;
	}

	public double getdevslope(){
		CurveFitter cf = new CurveFitter(x,y);
		cf.doFit(fitType);
		double[] aux = cf.getParams(); 							
		slope = aux[1];
		origin = aux[0];
		double Sx = 0;
		//double Sy = 0;
		//double Sxy = 0;
		double Sxx = 0;
		double chi=0;
		for (int i = 0; i <x.length; i++){
			Sx +=x[i];
			//Sy +=y[i];
			//Sxy +=x[i]*y[i];
			Sxx +=pow(x[i],2);
			chi += pow((y[i] - slope*x[i] - origin),2);
		}
		devslope = sqrt((x.length*chi)/((x.length*Sxx - Sx*Sx)*(x.length -2)));
		return devslope;
	}
	


}
				
