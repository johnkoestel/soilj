package soilj;

import ij.ImagePlus;
import ij.plugin.PlugIn;
import soilj.tools.ImageManipulator;
import soilj.tools.InputOutput;
import soilj.tools.MenuWaiter;

import java.io.File;

public class NormalizeTheseImages_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ImageManipulator jIM = new ImageManipulator();	
		MenuWaiter menu = new MenuWaiter();
		
		MenuWaiter.NormalizerReferences myNR;
		
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus nowTiff = new ImagePlus();		
		ImagePlus outTiff = new ImagePlus();
		
		//ask for threshold choice
		myNR = menu.showNormalizerReferencesMenu();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
	
		//read gauges in the gauge folder
		String myGaugeFolder = jIO.findTheInnerCircle(myBaseFolder);		
		String[] myGauges = jIO.listGaugeFiles(myGaugeFolder, "");
	
		String myOutFolder;
		String myOutPath = null;		
		myOutFolder = "Normalized2" + myNR.lowerTag + "And" + myNR.upperTag;		
		myOutPath = myBaseFolder + myOutFolder;
		new File(myOutPath).mkdir();		
				
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  //myTiffs.length

			//try to free up some memory
			System.gc();
			
			//load image
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffs[i]);			
			
			//select the correct gauge and surface files			
			int[] myGandS = jIO.getTheCorrectGaugeNSurfaceFiles(myTiffs[i], myGauges, null);	
					
			//apply segmentation			
			outTiff = jIM.normalizeImages2AirAndWallGrey(nowTiff, myGauges[myGandS[0]], myNR, myOutPath);			
		
			//save result
			jIO.tiffSaver(myOutPath, myTiffs[i], outTiff);
			
		}
		
		//finally also save the gauge path..
		jIO.writeStringIntoAsciiFile(myOutPath + "\\Path2Gauge.txt", myGaugeFolder);
		
	}
	
}