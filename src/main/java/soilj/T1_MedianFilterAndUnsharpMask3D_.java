package SoilJ;

import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import SoilJ.tools.ImageManipulator;
import SoilJ.tools.InputOutput;
import SoilJ.tools.MenuWaiter;

import java.io.File;

public class T1_MedianFilterAndUnsharpMask3D_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();		
				
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();
		ImagePlus outTiff = new ImagePlus();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		
		//open dialog box to query the filter settings
		MenuWaiter.MedianFilterAndUnsharpMaskReturn mMUS = menu.showMedianAndUsharpMaskMenu();
		
		//load list of images to process
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
		String myTiffName;
		String myOutFolder = "UnsharpMedian";
						
		//create output folder
		String mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder);
		String myOutPath = mySubBaseFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir();		
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  	//myTiffs.length

			//assign current TIFF
			myTiffName = myTiffs[i];
			
			//load TIFF and apply filters
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffName);						
			outTiff = jIM.applyMedianFilterAndUnsharpMask(nowTiff, mMUS);	
			
			//save result
			jIO.tiffSaver(myOutPath, myTiffName, outTiff);
			
			//try to free up some memory
			IJ.freeMemory();IJ.freeMemory();			
		}
	}
}