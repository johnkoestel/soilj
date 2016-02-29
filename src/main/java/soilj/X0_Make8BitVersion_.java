package SoilJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import SoilJ.tools.InputOutput;

import java.io.File;

public class X0_Make8BitVersion_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();		
	
		// init variables
		int i;
	
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));		
		 	
		//create output paths
		String myOutFolder = "8bit";		
		String myOutPath = myBaseFolder + myOutFolder;
		new File(myOutPath).mkdir();
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffs[i]);	
			ImagePlus outTiff = new ImagePlus();
			ImageStack outStack = new ImageStack(nowTiff.getWidth(), nowTiff.getHeight());
			
			for (int j = 1 ; j < nowTiff.getNSlices() + 1 ; j++) {
				
				nowTiff.setSlice(j);
				ImageProcessor nowIP = nowTiff.getProcessor();
				nowIP.max(1);
				nowIP.multiply(255);
				
				ImageProcessor binIP = nowIP.convertToByteProcessor(false);
				
				outStack.addSlice(binIP);
				
			}
			
			outTiff.setStack(outStack);
										
			//apply segmentation	
			jIO.tiffSaver(myOutPath, myTiffs[i], outTiff);
			
		}
	}
}