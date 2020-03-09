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
import java.io.File;
import static java.lang.Math.*;	

public class _2D_Lac_Conf_Entropy implements PlugInFilter {
	public String base_path;
	public String imagename;	
	public ImagePlus salida;
	public int choose;
	public int setup(String arg, ImagePlus imp) {
		// INFORMAR SOBRE COMO SE REALIZAN LOS CALCULOS. DAR LA OPCION DE INVERTIR LA IMAGEN O CANCELAR
		GenericDialog gd1	= new GenericDialog("Process image");			// Create dialog
		String[] titleArray2 	= new String[2];
		gd1.addMessage("Image must be binary. Mass will be the number of white boxes. Do you want to analyze this image or the inverted one?");
		titleArray2[0] = "This image";
		titleArray2[1] = "Inverted image";
		gd1.addChoice("  ", titleArray2, titleArray2[0]);
		gd1.showDialog();
		//if (gd1.wasCanceled()){return;}
		choose = gd1.getNextChoiceIndex();
		if (choose == 0){}
		else if (choose == 1){// Then invert image
			IJ.run("Invert LUT");
		}
		base_path = imp.getOriginalFileInfo().directory;
		imagename = imp.getOriginalFileInfo().fileName;


		return DOES_8G+DOES_16+DOES_32;	// Suports gray scale image
		
	}


	public void run(ImageProcessor ip) {
		//Input dialog por saving path
		String save_path_pred;
		String[] parts = imagename.split("\\.");
		String part1 = parts[0]; // 004
		base_path = base_path + part1 + "/";
		if (choose == 0){save_path_pred = base_path + "gliding/white_max_value/min_1px/";} 
		else {save_path_pred = base_path + "gliding/binary/black_max_value/min_1px/";}		

		//GenericDialog gdpath = new GenericDialog("Select path for saving results");	// Create dialog
		//gdpath.addStringField("Ruta", save_path_pred);
		//gdpath.showDialog();
		String path;
		//path = gdpath.getNextString();
		path = save_path_pred;

		// Check if path exists and if not, create it
		File wdir = new File(path);
		if (!wdir.isDirectory() || !wdir.exists() || !wdir.canRead()) {
			wdir.mkdirs();
			IJ.log("Path created.");}
		
		// Plot tables and result
		lac_conf_entropy_function lac_entropy_fun = new lac_conf_entropy_function();
		//String str = new String(" 2D Multifractal Box Counting method");
		//lac_entropy_fun.str = "hola";
		lac_entropy_fun.ip = ip;
		lac_entropy_fun.path = path;
		lac_entropy_fun.out();
	}
}
