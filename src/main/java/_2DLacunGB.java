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
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.io.File;

/**
 * Tool that computes Lacunarity analysis using Gliding Box sliding method
 * on 2D Gray images
 * 
 * @author Ivan G Torre
 */

public class _2DLacunGB implements PlugInFilter {	
	public String base_path;
	public String imagename;	
	public ImagePlus salida;
	public int choose;
	
	/** ----------Choose to invert the image------------- */         
	public int setup(String arg, ImagePlus imp) {
		IJ.log(" 	MULTIFRAC  Copyright (C) <2020>  <ivangtorre>\n" +
				"	License GNU General Public License 3\n"+
				"	This program comes with ABSOLUTELY NO WARRANTY\n" + 
				"   This is free software, and you are welcome to redistribute it\n" + 
				"   under certain conditions. Please cite:\n" +
				"	I. G. Torre and A. M. Tarquis, MULTIFRAC: An ImageJ plugin for multifractal and \n" + 
				"	multiscale characterization of 2D and 3D stack images");
		
		GenericDialog gd1	= new GenericDialog("Process image");
		String[] titleArray2 	= new String[2];
		gd1.addMessage("Maximum measure will be for the white. Do you want to analyze this image or the inverted one?");
		titleArray2[0] = "This image";
		titleArray2[1] = "Inverted image";
		gd1.addChoice("  ", titleArray2, titleArray2[0]);
		gd1.showDialog();
		choose = gd1.getNextChoiceIndex();
		if (choose == 0){}
		else if (choose == 1){IJ.run("Invert");}
		base_path = imp.getOriginalFileInfo().directory;
		imagename = imp.getOriginalFileInfo().fileName;
		
		return DOES_8G+DOES_16+DOES_32;
		
	}


	public void run(ImageProcessor ip) {
    	/** ---------Handle output path---------------------------------- */         
		String save_path_pred;
		String[] parts = imagename.split("\\.");
		String part1 = parts[0]; // 004
		base_path = base_path + part1 + "/";
		if (choose == 0){save_path_pred = base_path + "gliding/white_max_value/min_1px/";} 
		else {save_path_pred = base_path + "gliding/black_max_value/min_1px/";}

		String path;
		path = save_path_pred;

    	/** ---------Check if path exists and create if necessary---------- */         
		File wdir = new File(path);
		if (!wdir.isDirectory() || !wdir.exists() || !wdir.canRead()) {
			wdir.mkdirs();
			IJ.log("Path created.");}
		
    	/** ---------Run the analysis---------------------------------- */         
		ToolLacun lacunarity = new ToolLacun();
		lacunarity.ip = ip;
		lacunarity.path = path;
		lacunarity.out();
	}
}
