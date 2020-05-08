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
import ij.IJ;
import ij.gui.*;

/**
 * Handles dialog for which color is max (black or white)
 * 
 * @author Ivan G Torre
 */

public class ColorDialog{
	int colorchoose;
	public int askmax(){
		GenericDialog gd1	= new GenericDialog("Process image");
		String[] titleArray2 	= new String[2];
		gd1.addMessage("Occupied boxes wil be black pixels. Do you want to analyze this image or to invert it?");
		titleArray2[0] = "This image";
		titleArray2[1] = "Invert image";
		gd1.addChoice("  ", titleArray2, titleArray2[0]);
		gd1.showDialog();
		if (gd1.wasCanceled()) return colorchoose;
		colorchoose = gd1.getNextChoiceIndex();
		if (colorchoose == 0){
		}
		else if (colorchoose == 1){
			IJ.run("Invert");
			IJ.wait(500);
		}
		return colorchoose;
	}
}
