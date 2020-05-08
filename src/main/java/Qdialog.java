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
import ij.gui.*;

/**
 * Handles dialog for q range input
 * 
 * @author Ivan G Torre
 */

public class Qdialog{
	double qmin, qmax, qincrement;
	
	public double[] getq(){
		
		GenericDialog gd = new GenericDialog("Select q values for multifractal dimensions");
		gd.addNumericField("qmin: ", -10, 1);
		gd.addNumericField("qmax: ", 10, 0);
		gd.addNumericField("qincrement: ", 1, 1);
		gd.showDialog();
		qmin       = (double)gd.getNextNumber();
		qmax       = (double)gd.getNextNumber();
		qincrement = (double)gd.getNextNumber();
		double[] q = new double[(int)((qmax - qmin)/qincrement)+1];
		for (int count = 0; count<q.length; count++){
			q[count] = qmin + count*qincrement;
		}
		return q;
	}
}

