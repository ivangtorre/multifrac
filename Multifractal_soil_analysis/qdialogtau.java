// Import libraries
import ij.plugin.filter.PlugInFilter;
import ij.gui.*;
import static java.lang.Math.*;	


public class qdialogtau{
	double qmin, qmax, qincrement;
	public double[] getq(){
		
			GenericDialog gd = new GenericDialog("Select q values for multifractal dimensions");	// Create dialog
			gd.addNumericField("qmin: ", 0.1, 1);							// Field1
			gd.addNumericField("qmax: ", 5, 0);							// Field2
			gd.addNumericField("qincrement: ", 0.1, 1);						// Field3
			gd.showDialog();									// Show
			qmin       = (double)gd.getNextNumber();						// Save value1
			qmax       = (double)gd.getNextNumber();						// Save value2
			qincrement = (double)gd.getNextNumber();						// Save value3
			//if (gd.wasCanceled())									// If no input, cancell program
	
		//while ( ((qmax - qmin) % qincrement) != 0);							// If not properly q values, repeat input

		double[] q = new double[(int)((qmax - qmin)/qincrement)+1];					// Declare q array
		for (int count = 0; count<q.length; count++){							// Asign q values to the array	
			q[count] = qmin + count*qincrement;
		}
		return q;
	}
}

