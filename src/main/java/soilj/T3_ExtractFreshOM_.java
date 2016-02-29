package soilj;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import soilj.tools.InputOutput;
import soilj.tools.MenuWaiter;
import soilj.tools.ObjectDetector;

import java.io.File;

public class T3_ExtractFreshOM_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		MenuWaiter menu = new MenuWaiter();
		
		// init variables
		int i;
		MenuWaiter.OMFinderSettings oMF = menu.new OMFinderSettings();
		oMF = menu.showOMFinderMenu();
		
		//construct image related objects
		ImagePlus nowTiff;  //input
		ImagePlus outTiff;  //output
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
		String myTiffName;
		String myOutFolder = "FreshOrganicMaterial";
						
		//create output folder
		String mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder);
		String myOutPath = mySubBaseFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir();
				
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  

			//assign current TIFF
			myTiffName = myTiffs[i];
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffName);
			
			//do the cutting
			outTiff = jOD.extractFreshOM(nowTiff, oMF);
			
			//save the results
			jIO.tiffSaver(myOutPath, myTiffName, outTiff);			
			
		}
	}
}