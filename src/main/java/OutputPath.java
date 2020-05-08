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
import java.io.File;

import ij.IJ;

/**
 * Handles output path
 * 
 * @author Ivan G Torre
 */

public class OutputPath{
String imagename;
String base_path;
int colorchoose;
	
	public String savepath(){
		String save_path_pred;
		IJ.log(""+imagename);
		String[] parts = imagename.split("\\.");
		String part1 = parts[0]; // 004
		base_path = base_path + part1 + "/";
		if (colorchoose == 0){save_path_pred = base_path + "box_counting/black_max_value/";} 
		else {save_path_pred = base_path + "box_counting/white_max_value/";}
		//String path;
		//path = save_path_pred;
		
		/** ---------Check if path exists and if not, create it----------------------------- */         
		File wdir = new File(save_path_pred);
		if (!wdir.isDirectory() || !wdir.exists() || !wdir.canRead()) {
		    wdir.mkdirs();
		    IJ.log("Path created.");}

		return save_path_pred;
	}
}


