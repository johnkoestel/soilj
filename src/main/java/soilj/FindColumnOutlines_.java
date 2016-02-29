package soilj;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.plugin.PlugIn;
import soilj.tools.ImageManipulator;
import soilj.tools.InputOutput;
import soilj.tools.MenuWaiter;
import soilj.tools.ObjectDetector;

import java.io.File;

public class FindColumnOutlines_ extends ImagePlus implements PlugIn  {

	// loads all 3D TIFFs in the specified folder and puts them straight
	// and aligns them as such that they stand 'perfectly' straight

	public void run(String arg) {

		//fix memory settings		
		Prefs.keepUndoBuffers=false;		
		
		//construct biggish objects
		InputOutput jIO = new InputOutput();		
		ObjectDetector jOD = new ObjectDetector();
		ImageManipulator jIM = new ImageManipulator();
		MenuWaiter menu = new MenuWaiter();
		
		ObjectDetector.ColumnContainer colCon = jOD.new ColumnContainer();
				
		// init variables
		int i;
		
		//construct image related objects
		ImagePlus nowTiff;  //perfectly upright columns
					
		//get details on PVC column finding
		MenuWaiter.ColumnFinderMenuReturn jCFS = menu.new ColumnFinderMenuReturn(); 
		jCFS = menu.showColumnFinderDialog();
		
		//read base folder and number of 3D images
		String myBaseFolder = jIO.chooseAFolder("Please choose the folder with your image data");
		String[] myTiffs = jIO.listTiffsInFolder(new File(myBaseFolder));
		String myTiffName;
		String myOutFolder = "WallCoordinateIdentified";
						
		//create output folder
		String mySubBaseFolder = jIO.getTheFolderAbove(myBaseFolder);
		String myOutPath = mySubBaseFolder + "\\" + myOutFolder;
		new File(myOutPath).mkdir();
		
		//create output folder for the binary mask of just the soil
		String myGaugeFolder = "InnerCircle";
		String myGaugePath = myOutPath + "\\" + myGaugeFolder;
		new File(myGaugePath).mkdir();
		
		//loop over 3D images
		for (i = 0 ; i < myTiffs.length ; i++) {  

			//assign current TIFF
			myTiffName = myTiffs[i];
			nowTiff = jIO.openTiff3D(myBaseFolder + "\\" + myTiffName);
			
		
			//do the cutting
			colCon = jIM.innerCircleAndTrimAluOrPVCColumns(nowTiff, jCFS);
				
			//save the results
			jIO.tiffSaver(myOutPath, myTiffName, colCon.cutTiff);
			jIO.writeEllipsoidMaskFile(myGaugePath + "\\Gauge_" + myTiffName.subSequence(0, myTiffName.length() - 4) + ".txt", colCon.mendedCC);
				
			//free memory
			colCon.nowTiff.flush();colCon.cutTiff.flush();
			
			nowTiff.flush();			
			IJ.freeMemory();IJ.freeMemory();
			
		}
	}
}