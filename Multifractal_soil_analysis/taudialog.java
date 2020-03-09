// Import libraries
import ij.plugin.filter.PlugInFilter;
import ij.gui.*;
import static java.lang.Math.*;	


public class taudialog{
	int taumin, taumax, tauincrement;
	public int[] gettau(){
		
		// DIALOG
		GenericDialog gd = new GenericDialog("Select Tau values for getting slope");	// Create dialog
		gd.addNumericField("taumin: ", 1, 0);						// Field1
		gd.addNumericField("taumax: ", 10, 0);						// Field2
		gd.addNumericField("tauincrement: ", 1, 0);					// Field3
		gd.showDialog();								// Show
		taumin       = (int)gd.getNextNumber();						// Save value1
		taumax       = (int)gd.getNextNumber();						// Save value2
		tauincrement = (int)gd.getNextNumber();						// Save value3

		// CREATE AND RETURN LIST
		int[] taulist = new int[(int)((taumax - taumin)/tauincrement)+1];		// Tau list
		for (int count = 0; count<taulist.length; count++){				// Asign Tau to list	
			taulist[count] = taumin + count*tauincrement;
		}
		return taulist;
	}
}

